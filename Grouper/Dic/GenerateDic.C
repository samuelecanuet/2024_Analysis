
#include <vector>
#include <utility> 


void GenerateDic()
{
    gSystem->AddDynamicPath("/home/local1/Documents/lib/SignalDict/");
    gSystem->Load("SignalDict.so");
    gSystem->AddIncludePath("-I/home/local1/Documents/lib/SignalDict/");
    // gSystem->AddIncludePath("home/local1/Documents/lib/SignalDict/Signal.h");
    // gSystem->AddIncludePath("../../../lib/SignalDict/");
    gInterpreter->GenerateDictionary("vector<vector<std::pair<Signal, Signal>>>", "vector;utility;Signal.h");
}