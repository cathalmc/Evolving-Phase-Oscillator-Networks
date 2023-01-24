#include <iostream>
#include <execution>
#include <vector>
#include <fstream>
#include <random>
#include <stdio.h>
#include <string>
#include <sstream>
#include <chrono> 



constexpr int nIterations;
constexpr int nInput = 5;
constexpr int inputNeurons[nInput] = { 0,1,2,3,4 };
constexpr int leftNeurons[nInput] = { 10,11,12,13,14 };
constexpr int rightNeurons[nInput] = { 20,21,22,23,24 };

const std::string filePath;

// = "C:\\Users\\catha\\Documents\\Research\\Models\\Other experiments\\AGI\\";

std::random_device dev;
std::mt19937 rng(dev());

#define LOG(x)  std::cout << x << std::endl //short hand for print function

template <typename T>
void writeBinaryFile(T* a, int length, std::string fileName, std::string location) {
    //writes a binary file to the output folder in "location"
    //a is a pointer to the memory where the data is located, length is the number of items
    //binary files are faster to read and write, and can be read in python with the numpy function "fromfile"
    //std::string folder = "popData\\";
    std::string fileType = ".bin";
    //location.append(folder);
    location.append(fileName);
    location.append(fileType);

    std::ofstream myfile;
    myfile.open(location, std::ios::out | std::ios::binary);

    if (myfile.is_open()) {
        myfile.write(reinterpret_cast<char*>(a), sizeof(T) * length);
        //LOG("File saved: " << fileName);
    }
    else {
        LOG("Could not open file: " << fileName);
        LOG("Full " << location);
    }
}

template <typename T>
T* readBinaryFile(std::string fileName, std::string location) {
    //reads binary file from output folder of "location"
    //std::string folder = "popData\\";
    std::string fileType = ".bin";
    //location.append(folder);
    location.append(fileName);
    location.append(fileType);

    std::ifstream is(location, std::ios::in | std::ios::binary);
    if (is) {
        // get length of file:
        is.seekg(0, is.end);
        int length = is.tellg();
        is.seekg(0, is.beg);
        char* buffer = new char[length];
        is.read(buffer, length);
        if (!is) std::cout << "error: only " << is.gcount() << " could be read" << std::endl;
        T* dat = reinterpret_cast<T*>(buffer);
        is.close();
        //LOG("Read file " << location);
        return dat;
    }
    else {
        LOG("Could not read file: " << location);
        return NULL;
    }
}

template <typename T>
class dataLogger {
    // helper class for saving results
public:
    unsigned int nitems;
    unsigned int nsamples;
    std::string name;
    std::string filePath;
    T* dat;
    bool alreadySaved;
    dataLogger(unsigned int nitems, unsigned int nsamples, std::string name, std::string filePath) :
        nitems(nitems),
        nsamples(nsamples),
        name(name),
        filePath(filePath),
        alreadySaved(false)
    {
        dat = (T*)calloc(nitems * nsamples, sizeof(T));
    }
    void add(unsigned int itemInd, unsigned int timeInd, T data) { dat[timeInd * nitems + itemInd] = data; }
    void save() {
        if (alreadySaved) {
            LOG("trying to save previously saved net");
            std::abort();
        }
        else {
            writeBinaryFile(dat, nitems * nsamples, name, filePath);
            free(dat);
            alreadySaved = true;
        }
    }

};



class lPRC {
public:
    int x, op, tau, period, ip;
    lPRC(int P, int tau) :x(0), op(0), tau(tau), period(P), ip(0) {}
    void reset() {
        x = 0; op = 0; ip = 0;
    }
    void update() {
        op = 0;
        x+=  ip +1;
        if (x >= period) { //spike and reset
            op = tau;
            x = 0;
        }

        if (x < 0) x = 0;

        ip = 0;
    }

};

class Edge {
public:
    int from, to, w;
    Edge(int from, int to, int w) : from(from), to(to), w(w) {};
};

inline float hillToClimb(float x) {
    return 100 - (x - 10) * (x - 10);
}

inline int inputScaler(float ip) {
    float x = 0.05f * (ip - 50.0f);
    float fastsigmoid = 0.5 + 0.5 * x / (1.0f + abs(x));
    return (int)5 * (1 - fastsigmoid);
}

template <typename T>
T convertInput(char ip[]) {
    //converts command line arguments to the appropriate type (float, int etc) using templates
    std::istringstream ss(ip);
    T x;
    if (!(ss >> x)) {
        std::cerr << "Invalid input: " << ip << '\n';
    }
    else if (!ss.eof()) {
        std::cerr << "Trailing characters after input: " << ip << '\n';
    }
    return x;
}


void runSim(int simNo, std::vector<float> &costData);

int main(int argc, char** argv) {
    int nSims = convertInput<int>(argv[1]);
    std::vector<float> costData(nSims,-1.0f);
    std::vector<int> nums(nSims);
    for (int i = 0; i < nSims; i++) nums[i] = i;

    auto startTime = std::chrono::high_resolution_clock::now();

    std::for_each(std::execution::par,  nums.begin(), nums.end(), [&](int i) {
         runSim(i, costData);
        });
    
    auto endTime = std::chrono::high_resolution_clock::now();


    double timeTaken =
        std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime).count();
    timeTaken *= 1e-9;


    //LOG("Time taken: " << timeTaken);
    std::cout.precision(14);
    std::cout << timeTaken << ",";
    for (auto& v : costData) std::cout << (float)v << ",";
        
    return 0;
    //poss.save();
    //costs.save();
    //nodeActs.save();
    //LOG("done");
    //system("PAUSE");
}

void runSim(int simNo, std::vector<float> &costData) {
    std::string simName = "simData";
    simName.append(std::to_string(simNo));
    //LOG("Running sim " << simNo);
    int* simData = readBinaryFile<int>(simName, filePath);

    int nNodes = simData[0];
    int nEdges = simData[1];
    int nodeOffset = 2;
    int edgeOffset = nodeOffset + 2 * nNodes;
    std::vector<std::vector<int>> nodes(nNodes, std::vector<int>(2));
    std::vector<std::vector<int>> edges(nEdges, std::vector<int>(3));

    std::uniform_int_distribution<> randIC(0, 10000);


    for (int i = 0; i < nNodes; i++) {
        nodes[i][0] = simData[nodeOffset + 2 * i];
        nodes[i][1] = simData[nodeOffset + 2 * i + 1];
    }
    for (int i = 0; i < nEdges; i++) {
        edges[i][0] = simData[edgeOffset + 3 * i];
        edges[i][1] = simData[edgeOffset + 3 * i + 1];
        edges[i][2] = simData[edgeOffset + 3 * i + 2];
    }

    //for (auto n : nodes) LOG(n[0] << " " << n[1]);

    //vars
    float TS = 0.1f;
    float cost = 0.0f;
    

    //dataLogger<float> poss(1, nIterations, "poss", filePath);
    //dataLogger<float> costs(1, nIterations, "cost", filePath);
    //dataLogger<float> nodeActs(10, nIterations, "nodeActs", filePath);

    std::vector<lPRC> neurons;
    std::vector<Edge> synapses;
    for (auto& n : nodes) {
        neurons.push_back(lPRC(n[0], n[1]));
    }
    for (auto& e : edges) {
        synapses.push_back(Edge(e[0], e[1], e[2]));
    }

    float left = 0.0f, right = 0.0f;

    float position;

    std::vector<float> initialPositions;
    for (int i = -15; i < 35; i += 4) initialPositions.push_back((float)i); 
    float nInitialPos = (float)initialPositions.size();

    for (float& iP : initialPositions) {
        position = iP;
        
        for (auto& neur : neurons) {
            neur.reset();
            neur.x = randIC(dev) % neur.period;
        }
        for (int it = 0; it < nIterations; it++) {
            for (auto& neur : neurons) neur.update();

            for (int i = 0; i < nEdges; i++) {
                auto s = synapses[i];
                neurons[s.to].ip += s.w * (neurons[s.from].op);
            }

            left = 0.0f; right = 0.0f;
            for (int i = 0; i < nInput; i++) {
                neurons[inputNeurons[i]].ip += inputScaler(hillToClimb(position));
                left += (float)((float)neurons[leftNeurons[i]].op) / ((float)nInput);
                right += (float)((float)neurons[rightNeurons[i]].op) / ((float)nInput);
            }
            position += TS * (right - left);
            cost += ((position - 10.0f) * (position - 10.0f))/((float)(nIterations)) /nInitialPos; 
            
        }
    }
    free(simData);
    costData[simNo] = cost;///((float)nIterations);
    //std::cout << cost << ",";

}

/*
for (auto& initPos : initialPoss) {
        position = initPos;
        for (auto& neur : neurons) neur.reset();

        for (int it = 0; it < nIterations; it++) {
            for (auto& neur : neurons) neur.update();

            for (int i = 0; i < nEdges; i++) {
                auto s = synapses[i];
                neurons[s.to].ip += s.w * (neurons[s.from].op);
            }

            left = 0.0f; right = 0.0f;
            for (int i = 0; i < nInput; i++) {
                neurons[inputNeurons[i]].ip += inputScaler(hillToClimb(position));
                left += neurons[leftNeurons[i]].op / ((float)nInput);
                right += neurons[rightNeurons[i]].op / ((float)nInput);
            }

            position += TS * (right - left);

            cost += (((float)it / (float)nIterations) * ((float)it / (float)nIterations)) * (abs(position - 10.0f));
        }
    }

*/

    /*
    dataLogger<int> neurX(1,nIterations,"neurX",filePath);
    dataLogger<int> neurIp(1, nIterations, "neurip", filePath);
    dataLogger<int> neurOp(1, nIterations, "neurop", filePath);
    
    lPRC neur(20, 2);

    LOG(neur.period);
    LOG(neur.tau);
    
    for (int i = 0; i < nIterations; i++) {
        neur.ip = 1;
        neur.update();
        neurX.add(0, i, neur.x);
        neurIp.add(0, i, neur.ip);
        neurOp.add(0, i, neur.op);

    }

    
    neurX.save();
    neurIp.save();
    neurOp.save();
    */


/*

    int* simData = readBinaryFile<int>("simData0", filePath);
    int nNodes = simData[0];
    int nEdges = simData[1];
    int nodeOffset = 2;
    int edgeOffset = nodeOffset + 2 * nNodes;
    std::vector<std::vector<int>> nodes(nNodes, std::vector<int>(2));
    std::vector<std::vector<int>> edges(nEdges, std::vector<int>(3));

    std::uniform_int_distribution<> randIC(0, 100);


    for (int i = 0; i < nNodes; i++) {
        nodes[i][0] = simData[nodeOffset + 2 * i];
        nodes[i][1] = simData[nodeOffset + 2 * i + 1];
    }
    for (int i = 0; i < nEdges; i++) {
        edges[i][0] = simData[edgeOffset + 3 * i];
        edges[i][1] = simData[edgeOffset + 3 * i + 1];
        edges[i][2] = simData[edgeOffset + 3 * i + 2];
    }

    //for (auto n : nodes) LOG(n[0] << " " << n[1]);

    //vars
    float TS = 0.1f;
    float cost = 0.0f;
    int nIterations = 1000;

    dataLogger<float> poss(1, nIterations, "poss", filePath);
    dataLogger<float> costs(1, nIterations, "cost", filePath);
    dataLogger<float> nodeActs(10, nIterations, "nodeActs", filePath);

    std::vector<lPRC> neurons;
    std::vector<Edge> synapses;
    for (auto& n : nodes) {
        neurons.push_back(lPRC(n[0], n[1]));
    }
    for (auto& e : edges) {
        synapses.push_back(Edge(e[0], e[1], e[2]));
    }

    constexpr int nInput = 5;
    constexpr int inputNeurons[nInput] = { 0,1,2,3,4 };
    constexpr int leftNeurons[nInput] = { 10,11,12,13,14 };
    constexpr int rightNeurons[nInput] = { 20,21,22,23,24 };

    float left = 0.0f, right = 0.0f;

    float position;

    std::vector<float> initialPositions;
    for (int i = -15; i < 35; i += 4) initialPositions.push_back((float)i); //25


    for (auto& iP : initialPositions) {
        position = iP;
        for (auto& neur : neurons) {
            neur.reset();
            neur.x = 1;//randIC(dev) % neur.period;
        }
        for (int it = 0; it < nIterations; it++) {
            for (auto& neur : neurons) neur.update();

            for (int i = 0; i < nEdges; i++) {
                auto s = synapses[i];
                neurons[s.to].ip += s.w * (neurons[s.from].op);
            }

            left = 0.0f; right = 0.0f;
            for (int i = 0; i < nInput; i++) {
                neurons[inputNeurons[i]].ip += inputScaler(hillToClimb(position));
                left += neurons[leftNeurons[i]].op / ((float)nInput);
                right += neurons[rightNeurons[i]].op / ((float)nInput);
            }
            position += TS * (right - left);
            cost += (((float)it / (float)nIterations) * ((float)it / (float)nIterations)) * (abs(position - 10.0f));

            //poss.add(0, it, position);
            //costs.add(0, it, cost);
            //for (int i = 0; i < 10; i++) nodeActs.add(i, it, neurons[5+i].op);

        }
    }

    std::cout << cost << ",";
    std::cout << cost + 1.0f;
*/