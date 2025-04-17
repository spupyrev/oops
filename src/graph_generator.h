#pragma once

#include <string>
#include <vector>

using EdgeTy = std::pair<int, int>;

void genByClass(const std::string& graphClass, int& n, std::vector<EdgeTy>& edges);

void genComplete(std::vector<EdgeTy>& edges, int n);
void genCompleteBipartite(std::vector<EdgeTy>& edges, int n, int m);
void genCirculant(std::vector<EdgeTy>& edges, int n, const std::vector<int>& S);
