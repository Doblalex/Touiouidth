#pragma once

#include "headers.hpp"

struct PartitionElement{
    unordered_set<int> elements;
    unordered_set<int> refined;
    PartitionElement* next = NULL;
    PartitionElement* last = NULL;
};

struct PartitionRefinement{
    unordered_map<int, PartitionElement*> at;
    PartitionElement* start;
    PartitionElement* end;

    PartitionRefinement(int n);
    PartitionRefinement(vector<int> V);

    void refine(vector<int> arr);
    void print();
};