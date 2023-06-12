#include "partition.hpp"

PartitionRefinement::PartitionRefinement(int n)
{
    this->start = new PartitionElement();
    this->end = this->start;
    for (int i = 0; i < n; i++)
    {
        this->start->elements.insert(i);
        at[i] = this->start;
    }
}

PartitionRefinement::PartitionRefinement(vector<int> V)
{
    int n = V.size();
    this->start = new PartitionElement();
    this->end = this->start;
    for (int i = 0; i < n; i++)
    {
        this->start->elements.insert(V[i]);
        at[V[i]] = this->start;
    }
}

void PartitionRefinement::refine(vector<int> arr)
{
    vector<PartitionElement *> refinedEls;
    for (auto x : arr)
    {
        PartitionElement *el = at[x];
        if (el->refined.size() == 0)
            refinedEls.push_back(el);
        el->refined.insert(x);
        el->elements.erase(x);
    }
    for (auto el : refinedEls)
    {
        PartitionElement *newel = new PartitionElement();
        newel->elements = unordered_set<int>(el->refined.begin(), el->refined.end());
        el->refined.clear();
        newel->last = this->end;
        newel->last->next = newel;
        this->end = newel;
        for (auto x : newel->elements)
            at[x] = newel;

        if (el->elements.size() == 0)
        {
            if (el->next == NULL)
            {
                el->last->next = NULL;
                this->end = el->last;
            }
            else if (el->last == NULL)
            {
                el->next->last = NULL;
                this->start = el->next;
            }
            else
            {
                el->last->next = el->next;
                el->next->last = el->last;
            }
            delete el;
        }
    }
}

void PartitionRefinement::print()
{
    auto el = this->start;
    while (el != NULL)
    {
        for (auto x : el->elements)
        {
            cout << x << " ";
        }
        cout << endl;
        el = el->next;
    }
}