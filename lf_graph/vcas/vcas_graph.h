#ifndef GRAPH_H
#define GRAPH_H

#include "rq_provider.h"
#include "random.h"
#include "plaf.h"

namespace vcas_graph {
using namespace std;

/////////////////////////////////////////////////////////////////
// TYPES
/////////////////////////////////////////////////////////////////

// ENode structure
// template <typename N>
// class node_t {
//    N * node;
// }
// node_t<vertex_t<K,V>>

template <typename K, typename V>
class vertex_t {
public:
    union {
        struct {
            K key;
            V val;
	         edge_t<K>* ehead; // pointer to the EHead
            vcas_obj_t<vertex_t<K, V>*> *volatile vnext; // pointer to the next VNode
        };
    };
};

template <typename K, typename V>
class edge_t {
public:
    union {
        struct {
            K key;
            vertex_t<K, V>* pointv; // pointer to its vertex
	         vcas_obj_t<edge_t<K>*> *volatile enext; // pointer to the next ENode
        };
    };
};

#define vertex_ptr vertex_t<K, V>*;
#define edge_ptr edge_t<K>*;


template <typename K, typename V, class RecManager>
class graph {
 private:
    volatile char padding0[PREFETCH_SIZE_BYTES]; // TODO: what is this for
    vertex_ptr volatile *vlistHead;
    vertex_ptr volatile *vlistTail;
    volatile char padding2[PREFETCH_SIZE_BYTES];

    const int NUM_PROCESSES;
    RecManager* const recmgr;
    Random* const threadRNGs;  // threadRNGs[tid * PREFETCH_SIZE_WORDS] = rng for thread tid
    RQProvider<K, V, graph<K, V, RecManager>, RecManager, true, false>* rqProvider;
 #ifdef USE_DEBUGCOUNTERS
    debugCounters* const counters;
 #endif

    edge_ptr allocateEdge(const int tid);
    void initEdge(const int tid, edge_ptr p_node, K key);
    vertex_ptr allocateVertex(const int tid);
    void initVertex(const int tid, vertex_ptr p_node, K key, V value);

    // Find pred and curr for VNode(key)
    void locateVPlus(vertex_ptr startV, vertex_ptr* n1, vertex_ptr* n2, K key);

    // Contains++
    bool ContainsVPlus(vertex_ptr* n1, vertex_ptr* n2, K key1, K key2);

    // Find pred and curr for ENode(key)
    void locateEPlus(edge_ptr startE, edge_ptr* n1, edge_ptr* n2, K key);
    
    // Find pred and curr for VNode(key)     
    void locateCPlus(vertex_ptr startV, vertex_ptr* n1, vertex_ptr* n2, K key);

    // Contains++   
    bool ContainsCPlus(vertex_ptr* n1, vertex_ptr* n2, K key1, K key2);

 public:
    graph(const int numProcesses, const K _KEY_MIN, const K _KEY_MAX,
           const V NO_VALUE, Random* const threadRNGs);
    ~graph();

   // creation of new ENode
   edge_ptr AllocateE(K key);

   // creation of new VNode
   vertex_ptr AllocateV(K key, V val);
 
   // init of sentinals Head and Tail
   void init();

   // TODO: this is not initialized well I do not believe ...
   void initGraph(int n) {
      int i, j;
      for (i = 1; i <= n; i++) {
         AddV(i);
      }

      for (i = 1; i <= n; i += 2) {
         for (j = i + 1; j <= n; j += 2) {
            int u = rand() % n + 1;
            int v = rand() % n + 1;
            if (u != v)
               AddE(u,v);
         }
      }
   }
    
   // add a new vertex in the vertex-list
   bool AddV(K key, V val);
    
   // Deletes the vertex from the vertex-list
   bool RemoveV(K key);
   
   // Contains VNode
   bool ContainsV(K key);

   // add a new edge in the edge-list
   bool AddE(K key1, K key2); 
        
   // Deletes an edge from the edge-list if present
   bool RemoveE(K key1, K key2);
 
   //Contains ENode       
   bool ContainsE(K key1, K key2);
};
} // namespace
#endif  // GRAPH_H