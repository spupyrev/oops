#pragma once

#include "adjacency.h"
#include "io_graph.h"

#include <functional>
#include <string>
#include <vector>

struct InputGraph;
struct Result;

std::vector<std::pair<std::string, AdjListTy>> readMagma(const std::string& filename, const int maxN);
std::vector<std::pair<std::string, AdjListTy>> readLst(const std::string& filename, const int maxN);

std::vector<std::pair<std::string, AdjListTy>> readCfg(
    const std::string& filename, 
    const std::string& part,
    const std::function<bool(const int, const std::vector<EdgeTy>&)> graphFilter);

std::vector<std::pair<std::string, AdjListTy>> readG6(
    const std::string& filename, 
    const std::string& part,
    const std::function<bool(const int, const std::vector<EdgeTy>&)> graphFilter);

std::vector<std::pair<std::string, AdjListTy>> readS6(
    const std::string& filename, 
    const std::string& part,
    const std::function<bool(const int, const std::vector<EdgeTy>&)> graphFilter);


void printText(const std::string& filename, const std::string& graphName, const int n, const std::vector<EdgeTy>& edges);
void printGML(const std::string& filename, int n, const std::vector<EdgeTy>& edges);
void printInput(const std::string& filename, const std::string& graphName, const InputGraph& graph);
void printOutput(const std::string& filename, const std::string& graphName, const InputGraph& graph, const Result& result);
void printResultRaw(const int n, const std::vector<EdgeTy>& edges, const Result& result);