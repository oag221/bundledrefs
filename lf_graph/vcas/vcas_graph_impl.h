/*
 * File:lfGraphDS.cpp
 *
 * Author(s):
 *   Bapi Chatterjee <bapchatt@in.ibm.com>
 *   Sathya Peri <sathya_p@iith.ac.in>
 *   Muktikanta Sa   <cs15resch11012@iith.ac.in>
 *   Nandini Singhal <cs15mtech01004@iith.ac.in>
 *   
 * Description:
 *   lock-free implementation of a graph
 * Copyright (c) 2018.
 * last Updated: 31/07/2018
 *
*/
#ifndef GRAPH_IMPL_H
#define GRAPH_IMPL_H

#include <iostream>
#include <float.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <assert.h>
#include <getopt.h>
#include <limits.h>
#include <signal.h>
#include <sys/time.h>
#include <time.h>
#include <vector>
#include <ctime>        // std::time
#include <random>
#include <algorithm>
#include <iterator>
#include <math.h>
#include <time.h>
#include <fstream>
#include <iomanip>
#include <sys/time.h>
#include <atomic>
#include <list>
#include <queue>
#include <stack>

#include "vcas_graph.h"

// ofstream coutt("getpath.txt");
// freopen("error.txt", "w", stderr);
  
#define vntp "VERTEX NOT PRESENT"
#define entp "EDGE NOT PRESENT"
#define ventp "VERTEX OR EDGE NOT PRESENT"
#define ep "EDGE PRESENT"
#define eadd "EDGE ADDED"
#define er "EDGE REMOVED"
#define ef "EDGE FOUND"
#define vp "VERTEX PRESENT"
#define vadd "VERTEX ADDED"
#define vr "VERTEX REMOVED"

namespace vcas_graph {
using namespace std;


inline int is_marked_ref(long i) {
  return (int) (i & 0x1L);
}


inline long unset_mark(long i) {
  i &= ~0x1L;
  return i;
}

inline long set_mark(long i) {
  i |= 0x1L;
  return i;
}

inline long get_unmarked_ref(long w) {
  return w & ~0x1L;
}

inline long get_marked_ref(long w) {
  return w | 0x1L;
}

template <typename K, typename V, class RecManager>
graph<K, V, RecManager>::graph(const int numProcesses, const K _KEY_MIN,
                              const K _KEY_MAX, const V NO_VALUE,
                              Random* const threadRNGs)
   : NUM_PROCESSES(numProcesses),
      recmgr(new RecManager(numProcesses, 0)),
      threadRNGs(threadRNGs),
   #ifdef USE_DEBUGCOUNTERS
      counters(new debugCounters(numProcesses)),
   #endif
      KEY_MIN(_KEY_MIN),
      KEY_MAX(_KEY_MAX),
      NO_VALUE(NO_VALUE) {
   
   rqProvider = new RQProvider<K, V, graph<K, V, RecManager>, RecManager,
                     true, false>(numProcesses, this, recmgr);

   const int dummyTid = 0;
   recmgr->initThread(dummyTid);

   vlistHead = allocateVertex(dummyTid);
   vlistTail = allocateVertex(dummyTid);
   initVertex(dummyTid, vlistHead, _KEY_MIN, NO_VALUE);
   initVertex(dummyTid, vlistTail, _KEY_MAX, NO_VALUE);

   rqProvider->write_vcas(dummyTid, &vlistHead->vnext, vlistTail);
   rqProvider->write_vcas(dummyTid, &vlistTail->vnext, (vertex_ptr) nullptr);
}

// creation of new ENode
template <typename K, typename V, class RecManager>
edge_ptr graph<K, V, RecManager>::allocateEdge(const int tid) {
   edge_ptr newe = recmgr->template allocate<edge_t<K>>(tid); // TODO: does this exist ?? like w only one type in <K>
   if (newe == NULL) {
      cout << "ERROR: out of memory" << endl;
      exit(-1);
   }
   return newe;
}

// TODO: how to deal with marked being bits ?? -> I think it's fine as it is
template <typename K, typename V, class RecManager>
void graph<K, V, RecManager>::initEdge(const int tid, edge_ptr p_node, const K& key, vertex_ptr v_ptr) {
   rqProvider->init_node(tid, p_node);
   p_node->key = key;
   p_node->pointv = v_ptr;
   p_node->enext = new vcas_obj_t<edge_ptr>(nullptr, nullptr); // later we CAS to add this node to the structure
}

// creation of new VNode
template <typename K, typename V, class RecManager>
vertex_ptr graph<K, V, RecManager>::allocateVertex(const int tid) {
   vertex_ptr nnode = recmgr->template allocate<vertex_t<K, V>>(tid);
   if (nnode == NULL) {
    cout << "ERROR: out of memory" << endl;
    exit(-1);
  }
  return nnode;
}

template <typename K, typename V, class RecManager>
void graph<K, V, RecManager>::initVertex(const int tid, vertex_ptr newv, const K& key, const V& value) {
   rqProvider->init_node(tid, newv);
   // allocate and init sentinel nodes for edge list
   edge_ptr ehead = allocateEdge(tid);
   initEdge(tid, ehead, INT_MIN, newv);
   edge_ptr etail = allocateEdge(tid);
   initEdge(tid, etail, INT_MAX, newv);

   rqProvider->write_vcas(tid, &ehead->enext, etail);
   rqProvider->write_vcas(tid, &etail->enext, (edge_ptr) nullptr);

   newv->key = key;
   newv->val = value;
   newv->ehead = ehead;
   newv->vnext = new vcas_obj_t<edge_ptr>(nullptr, nullptr);
}

// Find pred and curr for VNode(key)
template <typename K, typename V, class RecManager>
void graph<K, V, RecManager>::locateVPlus(vertex_ptr startV, vertex_ptr* n1, vertex_ptr* n2, const K& key) {
   vertex_ptr succv, currv, predv;
   retry:
   while (true) {
      predv = startV;
      currv = rqProvider->read_vcas(tid, &predv->vnext);
      while (true) {
         succv = rqProvider->read_vcas(tid, &currv->vnext);
         while (rqProvider->read_vcas(tid, &currv->vnext) != NULL && is_marked_ref((long) succv) && currv->val < key ) {
            // TODO: do I need to make marked a field / vcas_obj or is this okay ??
            if (!rqProvider->cas_vcas(tid, &predv->vnext, currv, (vertex_ptr) get_unmarked_ref((long)succv)))
               goto retry;
            currv = (vertex_ptr) get_unmarked_ref((long)succv);
            succv = rqProvider->read_vcas(tid, &currv->vnext);
         }
         if (currv->val >= key) {
            (*n1) = predv;
            (*n2) = currv;
            return;
         }
         predv = currv;
         currv = succv;
      }
   }
}

// add a new vertex in the vertex-list
template <typename K, typename V, class RecManager>
pair<V,bool> graph<K, V, RecManager>::AddV(const int tid, const K& key, const V& val) { // TODO: what is "const K&""
   vertex_ptr predv, currv;
   pair<V,bool> ret = pair<V,bool>(NO_VALUE, false);
   recmgr->leaveQuiescentState(tid);
   while (true) {
      locateVPlus(vlistHead, &predv, &currv, key); // find the location, <pred, curr>
      if (currv->key == key) {
         ret = pair<V,bool>(currv->val, false);
         recmgr->enterQuiescentState(tid);
         return ret; // key already present
      } else {
         vertex_ptr newv = allocateVertex(tid); // create a new vertex node
         initVertex(newv);
         rqProvider->write_vcas(tid, &newv->vnext, currv);
         
         vertex_ptr insertedNodes[] = {newv, NULL};
         vertex_ptr deletedNodes[] = {NULL};
         rqProvider->linearize_update_at_cas(tid, &predv->vnext, currv, newv, insertedNodes, deletedNodes);
         ret = pair<V,bool>(NO_VALUE, true);
         break;
      }
   }
   recmgr->enterQuiescentState(tid);
   return ret;
}

// TODO: modify to include vcas in remaining methods
   
// Deletes the vertex from the vertex-list
template <typename K, typename V, class RecManager>
bool graph<K, V, RecManager>::RemoveV(const int tid, const K& key) {
   vertex_ptr predv, currv, succv;
   recmgr->leaveQuiescentState(tid);
   while (true) {
      locateVPlus(vlistHead, &predv, &currv, key);
      if (currv->val != key)
         recmgr->enterQuiescentState(tid);
         return false; // key is not present
      succv = rqProvider->read_vcas(tid, &currv->vnext);
      if (!is_marked_ref((long) succv)) {
         if (!rqProvider->cas_vcas(tid, &currv->vnext, succv, (vertex_ptr)get_marked_ref((long)succv))) // mark curr as logically deleted
            continue;
         if (rqProvider->cas_vcas(tid, &predv->vnext, currv, succv)) // change prev to point to the next vertex
            break;
      }
   }
   recmgr->enterQuiescentState(tid);
   return true;
}

// Contains VNode
template <typename K, typename V, class RecManager>
bool graph<K, V, RecManager>::ContainsV(const int tid, const K& key) {
   vertex_ptr currv = vlistHead;
   recmgr->leaveQuiescentState(tid);
   while (rqProvider->read_vcas(tid, &currv->vnext) && currv->val < key) {
      currv =  (vertex_ptr) get_unmarked_ref((long)rqProvider->read_vcas(tid, &currv->vnext));
   }
   vertex_ptr succv = rqProvider->read_vcas(tid, &currv->vnext);
   if (rqProvider->read_vcas(tid, &currv->vnext) && currv->val == key && !is_marked_ref((long) succv)) {
      recmgr->enterQuiescentState(tid);
      return true;
   } else {
      recmgr->enterQuiescentState(tid);
      return false;
   }   
}

// Contains++
template <typename K, typename V, class RecManager>
bool graph<K, V, RecManager>::ContainsVPlus(const int tid, vertex_ptr* n1, vertex_ptr* n2, const K& key1, const K& key2) {
   vertex_ptr curr1, pred1, curr2, pred2;
   if (key1 < key2) {
      locateVPlus(vlistHead, &pred1, &curr1, key1); //first look for key1
      if((!rqProvider->read_vcas(tid, &curr1->vnext)) || curr1->val != key1) {
         return false; // key1 is not present in the vertex-list
      }
      locateVPlus(curr1, &pred2, &curr2, key2); // looking for key2 only if key1 present
      if((!rqProvider->read_vcas(tid, &curr2->vnext)) || curr2->val != key2) {
         return false; // key2 is not present in the vertex-list
      }
         
   } else {
      locateVPlus(vlistHead, &pred2, &curr2, key2); //first look for key2 
      if((!rqProvider->read_vcas(tid, &curr2->vnext)) || curr2->val != key2) {
         return false; // key2 is not present in the vertex-list
      } 
      locateVPlus(curr2, &pred1, &curr1, key1); // looking for key1 only if key2 present
      if ((!rqProvider->read_vcas(tid, &curr1->vnext)) || curr1->val != key1) {
         return false; // key1 is not present in the vertex-list
      }   
   }
   (*n1) = curr1; 
   (*n2) = curr2; 
   return true;    
}

// Find pred and curr for ENode(key)
template <typename K, typename V, class RecManager>
void graph<K, V, RecManager>::locateEPlus(const int tid, edge_ptr startE, edge_ptr* n1, edge_ptr* n2, const K& key) {
   edge_ptr succe, curre, prede;
   vertex_ptr tv;
   
   retry: while (true) {
      prede = startE;

      curre = rqProvider->read_vcas(tid, &prede->enext);
      while (true) {
         succe = rqProvider->read_vcas(tid, &curre->enext);
         tv = rqProvider->read_vcas(tid, &curre->pointv);
         /*helping: delete one or more enodes whose vertex was marked*/
       retry2: while(tv && rqProvider->read_vcas(tid, &tv->vnext) && rqProvider->read_vcas(tid, &curre->enext) && is_marked_ref((long)rqProvider->read_vcas(tid, &tv->vnext))
                     && !is_marked_ref((long) succe) && curre->val < key) {
            if (!rqProvider->cas_vcas(tid, &curre->enext, succe, (edge_ptr)get_marked_ref((long)succe)))
               goto retry;
            if (!rqProvider->cas_vcas(tid, &prede->enext, curre, succe))
               goto retry;
            curre = (edge_ptr)get_unmarked_ref((long)succe);
            succe = rqProvider->read_vcas(tid, &curre->enext);
            tv = rqProvider->read_vcas(tid, &curre->pointv);
         }
         /*helping: delete one or more enodes which are marked*/
         while(rqProvider->read_vcas(tid, &curre->enext) && is_marked_ref((long) succe) && !is_marked_ref((long)rqProvider->read_vcas(tid, &tv->vnext)) &&  curre->val < key ){
            if(!rqProvider->cas_vcas(tid, &prede->enext, curre, (edge_ptr)get_unmarked_ref((long)succe)))
               goto retry;
            curre = (edge_ptr)get_unmarked_ref((long)succe);
            succe = rqProvider->read_vcas(tid, &curre->enext);
            tv = rqProvider->read_vcas(tid, &curre->pointv);
         }
      
         if(tv && rqProvider->read_vcas(tid, &tv->vnext) && is_marked_ref((long) rqProvider->read_vcas(tid, &tv->vnext)) && rqProvider->read_vcas(tid, &curre->enext) && curre->val < key)
            goto retry2;
         if(curre->val >= key){
            (*n1) = prede;
            (*n2) = curre;
            return;
         }
         prede = curre;
         curre =(edge_ptr)get_unmarked_ref((long)succe);
      }
   }
} 

// add a new edge in the edge-list
template <typename K, typename V, class RecManager>
bool graph<K, V, RecManager>::AddE(const int tid, const K& key1, const K& key2) {
   edge_ptr prede, curre;
   vertex_ptr u, v;
   recmgr->leaveQuiescentState(tid);
   bool flag = ContainsVPlus(&u, &v, key1, key2);
   if (flag == false)
      recmgr->enterQuiescentState(tid);
      return false; // either of the vertex is not present
   
   while (true) {
      if(is_marked_ref((long) rqProvider->read_vcas(tid, &u->vnext)) || is_marked_ref((long) rqProvider->read_vcas(tid, &v->vnext))) {
         recmgr->enterQuiescentState(tid);
         return false;
      }
         
      locateEPlus(rqProvider->read_vcas(tid, &u->enext), &prede, &curre, key2);
      if(curre->val == key2) {
         recmgr->enterQuiescentState(tid);
         return false; // edge already present 
      }
      
      edge_ptr newe = createE(key2);// create a new edge node

      rqProvider->write_vcas(tid, &newe->enext, curre); // connect newe->next to curr
      newe->pointv.store(v); // points to its vertex - note: store() bc immutable so not a vcas_obj

      if(rqProvider->cas_vcas(tid, &prede->enext, curre, newe)) { // actual insertion
         recmgr->enterQuiescentState(tid);
         return true;
      }
   }
}    
      
// Deletes an edge from the edge-list if present
template <typename K, typename V, class RecManager>
bool graph<K, V, RecManager>::RemoveE(const int tid, const K& key1, const K& key2) {
   edge_ptr prede, curre, succe;
   vertex_ptr u, v;
   recmgr->leaveQuiescentState(tid);
   bool flag = ContainsVPlus(&u, &v, key1, key2);
   if (flag == false) {
      recmgr->enterQuiescentState(tid);
      return false; // either of the vertex is not present
   }
   while (true) {
      if(is_marked_ref((long) rqProvider->read_vcas(tid, &u->vnext)) || is_marked_ref((long) rqProvider->read_vcas(tid, &v->vnext))) {
         recmgr->enterQuiescentState(tid);
         return false;
      }
      locateEPlus(rqProvider->read_vcas(tid, &u->enext), &prede, &curre, key2);
      if(curre->val != key2) {
         recmgr->enterQuiescentState(tid);
         return false; // edge already present
      }
      succe = rqProvider->read_vcas(tid, &curre->enext);
      if (!is_marked_ref((long) succe)) {
         if(!rqProvider->cas_vcas(tid, &curre->enext, succe, (edge_ptr)get_marked_ref((long)succe)))
            continue;
         if (!rqProvider->cas_vcas(tid, &prede->enext, curre, succe))
            break;
      }
   }
   // TODO: if we break out of while loop, what should be returned ?? --> is it possible right now that we mark a node as logically deleted but never physically delete it ?? (see two if statements a couple lines up)
   return false;
}
// Find pred and curr for VNode(key)
template <typename K, typename V, class RecManager>
void graph<K, V, RecManager>::locateCPlus(const int tid, vertex_ptr startV, vertex_ptr* n1, vertex_ptr* n2, const K& key) {
   vertex_ptr currv, predv;
   predv = startV;
   currv = rqProvider->read_vcas(tid, &startV->vnext);
   while (currv && currv->val < key) {
      predv = currv;
      currv = (vvertex_ptr)get_unmarked_ref((long)rqProvider->read_vcas(tid, &currv->vnext));
   }
   (*n1) = predv;
   (*n2) = currv;
   return;
}

// Contains++
template <typename K, typename V, class RecManager>
bool graph<K, V, RecManager>::ContainsCPlus(const int tid, vertex_ptr* n1, vertex_ptr* n2, const K& key1, const K& key2) {
   vertex_ptr curr1, pred1, curr2, pred2;
   if (key1 < key2) {
      locateCPlus(vlistHead, &pred1, &curr1, key1); //first look for key1
      if((!rqProvider->read_vcas(tid, &curr1->vnext)) || curr1->val != key1)
         return false; // key1 is not present in the vertex-list

      locateCPlus(curr1, &pred2, &curr2, key2); // looking for key2 only if key1 present
      if((!rqProvider->read_vcas(tid, &curr2->vnext)) || curr2->val != key2)
         return false; // key2 is not present in the vertex-list
   } else {
      locateCPlus(vlistHead, &pred2, &curr2, key2); //first look for key2 
      if((!rqProvider->read_vcas(tid, &curr2->vnext)) || curr2->val != key2)
         return false; // key2 is not present in the vertex-list
               
      locateCPlus(curr2, &pred1, &curr1, key1); // looking for key1 only if key2 present
      if((!rqProvider->read_vcas(tid, &curr1->vnext)) || curr1->val != key1)
         return false; // key1 is not present in the vertex-list
   }
   (*n1) = curr1;
   (*n2) = curr2;
   return true;
}

//Contains ENode
template <typename K, typename V, class RecManager>
bool graph<K, V, RecManager>::ContainsE(const int tid, const K& key1, const K& key2) {
   edge_ptr curre, prede;
   vertex_ptr u, v;
   recmgr->leaveQuiescentState(tid);
   bool flag = ContainsCPlus(&u, &v, key1, key2);
   if (flag == false) {
      recmgr->enterQuiescentState(tid);
      return false; // at least one of the vertices is not present
   }
   curre = rqProvider->read_vcas(tid, &u->enext);
   while(rqProvider->read_vcas(tid, &curre->enext) && curre->val < key2) {
      curre =  (edge_ptr)get_unmarked_ref((long)rqProvider->read_vcas(tid, &curre->enext));
   }
   if ((rqProvider->read_vcas(tid, &curre->enext)) && curre->val == key2 && !is_marked_ref((long)rqProvider->read_vcas(tid, &curre->enext)) && !is_marked_ref((long)rqProvider->read_vcas(tid, &u->vnext)) && !is_marked_ref((long) rqProvider->read_vcas(tid, &v->vnext))) {
      recmgr->enterQuiescentState(tid);
      return true;
   }
   recmgr->enterQuiescentState(tid);
   return false;
}

} // namespace
#endif /* GRAPH_IMPL_H */