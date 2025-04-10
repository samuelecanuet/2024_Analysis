
#include <vector>
#include <utility> 


void GenerateDic()
{
    gSystem->AddDynamicPath("/home/local1/Documents/lib/SignalDict/");
    gSystem->Load("SignalDict.so");
    gSystem->AddIncludePath("-I/home/local1/Documents/lib/SignalDict/");
    gInterpreter->GenerateDictionary("vector<vector<std::pair<Signal, Signal>>>", "vector;utility;Signal.h");

    gSystem->AddDynamicPath("/home/local1/Documents/lib/SignalDict2/");
    gSystem->Load("SignalDict2.so");
    gSystem->AddIncludePath("-I/home/local1/Documents/lib/SignalDict2/");
    gInterpreter->GenerateDictionary("vector<vector<std::pair<Signal2, Signal2>>>", "vector;utility;Signal2.h");
}