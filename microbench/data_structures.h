#ifndef DATA_STRUCTURE_H
#define DATA_STRUCTURE_H

const test_type NO_VALUE = -1;
const test_type KEY_MIN = numeric_limits<test_type>::min() + 1;
const test_type KEY_MAX =
    numeric_limits<test_type>::max() -
    1; // must be less than max(), because the snap collector needs a reserved
       // key larger than this!
#define KEY_PRECEEDING(key) (key - 1)

#ifdef RQ_SNAPCOLLECTOR
#define RQ_SNAPCOLLECTOR_OBJECT_TYPES                                      \
  , SnapCollector<node_t<test_type, test_type>, test_type>,                \
      SnapCollector<node_t<test_type, test_type>, test_type>::NodeWrapper, \
      ReportItem, CompactReportItem
#define RQ_SNAPCOLLECTOR_OBJ_SIZES                                             \
  << " SnapCollector="                                                         \
  << (sizeof(SnapCollector<node_t<test_type, test_type>, test_type>))          \
  << " NodeWrapper="                                                           \
  << (sizeof(                                                                  \
         SnapCollector<node_t<test_type, test_type>, test_type>::NodeWrapper)) \
  << " ReportItem=" << (sizeof(ReportItem))                                    \
  << " CompactReportItem=" << (sizeof(CompactReportItem))
#else
#define RQ_SNAPCOLLECTOR_OBJECT_TYPES
#define RQ_SNAPCOLLECTOR_OBJ_SIZES
#endif

#ifndef VCAS_GRAPH

#ifndef RQ_FUNC
#define RQ_FUNC rangeQuery
#endif

#ifndef INSERT_FUNC
#define INSERT_FUNC insert
#endif

#ifndef ERASE_FUNC
#define ERASE_FUNC erase
#endif

#ifndef FIND_FUNC
#define FIND_FUNC contains
#endif

#else
// graph specific functions
#ifndef INSERT_VERTEX_FUNC
#define INSERT_VERTEX_FUNC AddV
#endif

#ifndef REMOVE_VERTEX_FUNC
#define REMOVE_VERTEX_FUNC RemoveV
#endif

#ifndef FIND_VERTEX_FUNC
#define FIND_VERTEX_FUNC ContainsV
#endif

#ifndef INSERT_EDGE_FUNC
#define INSERT_EDGE_FUNC AddE
#endif

#ifndef REMOVE_EDGE_FUNC
#define REMOVE_EDGE_FUNC RemoveE
#endif

#ifndef FIND_EDGE_FUNC
#define FIND_EDGE_FUNC ContainsE
#endif
// end graph specific functions
#endif


#if defined ABTREE || defined BSLACK
#define VALUE ((void *)(int64_t)key)
#define KEY keys[0]
#define VALUE_TYPE void *
#else
#define VALUE key
#define KEY key
#define VALUE_TYPE test_type
#endif

#if defined ABTREE || defined BSLACK
#if defined ABTREE
#define USE_SIMPLIFIED_ABTREE_REBALANCING
#endif
#define ABTREE_DEGREE 16
#include "bslack_impl.h"
#include "record_manager.h"
using namespace bslack_ns;

#define DS_DECLARATION \
  bslack<ABTREE_DEGREE, test_type, less<test_type>, MEMMGMT_T>
#define MEMMGMT_T \
  record_manager<RECLAIM, ALLOC, POOL, Node<ABTREE_DEGREE, test_type>>
#define DS_CONSTRUCTOR \
  new DS_DECLARATION(TOTAL_THREADS, ABTREE_DEGREE, KEY_MAX, SIGQUIT)

// note: INSERT success checks use "== NO_VALUE" so that prefilling can tell
// that a new KEY has been inserted
#define INSERT_AND_CHECK_SUCCESS \
  ds->INSERT_FUNC(tid, key, VALUE) == ds->NO_VALUE
#define DELETE_AND_CHECK_SUCCESS ds->ERASE_FUNC(tid, key).second
#define FIND_AND_CHECK_SUCCESS ds->FIND_FUNC(tid, key)
#define RQ_AND_CHECK_SUCCESS(rqcnt)                               \
  (rqcnt) = ds->RQ_FUNC(tid, key, key + RQSIZE - 1, rqResultKeys, \
                        (VALUE_TYPE *)rqResultValues)
#define RQ_GARBAGE(rqcnt) rqResultKeys[0] + rqResultKeys[(rqcnt)-1]
#define INIT_THREAD(tid) ds->initThread(tid)
#define DEINIT_THREAD(tid) ds->deinitThread(tid)
#define INIT_ALL
#define DEINIT_ALL

#define PRINT_OBJ_SIZES                                                   \
  cout << "sizes: node=" << (sizeof(Node<ABTREE_DEGREE, test_type>))      \
       << " descriptor=" << (sizeof(SCXRecord<ABTREE_DEGREE, test_type>)) \
       << endl;

#elif defined(BST)
#include "bst_impl.h"
#include "record_manager.h"
using namespace bst_ns;

#define DS_DECLARATION bst<test_type, test_type, less<test_type>, MEMMGMT_T>
#define MEMMGMT_T \
  record_manager<RECLAIM, ALLOC, POOL, Node<test_type, test_type>>
#define DS_CONSTRUCTOR \
  new DS_DECLARATION(KEY_MAX, NO_VALUE, TOTAL_THREADS, SIGQUIT)

#define INSERT_AND_CHECK_SUCCESS \
  ds->INSERT_FUNC(tid, key, VALUE) == ds->NO_VALUE
#define DELETE_AND_CHECK_SUCCESS ds->ERASE_FUNC(tid, key).second
#define FIND_AND_CHECK_SUCCESS ds->FIND_FUNC(tid, key)
#define RQ_AND_CHECK_SUCCESS(rqcnt)                               \
  (rqcnt) = ds->RQ_FUNC(tid, key, key + RQSIZE - 1, rqResultKeys, \
                        (VALUE_TYPE *)rqResultValues)
#define RQ_GARBAGE(rqcnt) rqResultKeys[0] + rqResultKeys[(rqcnt)-1]
#define INIT_THREAD(tid) ds->initThread(tid)
#define DEINIT_THREAD(tid) ds->deinitThread(tid)
#define INIT_ALL
#define DEINIT_ALL

#define PRINT_OBJ_SIZES                                          \
  cout << "sizes: node=" << (sizeof(Node<test_type, test_type>)) \
       << " descriptor=" << (sizeof(SCXRecord<test_type, test_type>)) << endl;

#elif defined(CITRUS)
#include "record_manager.h"
#include "citrus_impl.h"

#define DS_DECLARATION citrustree<test_type, test_type, MEMMGMT_T>
#define MEMMGMT_T \
  record_manager<RECLAIM, ALLOC, POOL, node_t<test_type, test_type>>
#define DS_CONSTRUCTOR new DS_DECLARATION(MAXKEY, NO_VALUE, TOTAL_THREADS)

#define INSERT_AND_CHECK_SUCCESS \
  ds->INSERT_FUNC(tid, key, VALUE) == ds->NO_VALUE
#define DELETE_AND_CHECK_SUCCESS ds->ERASE_FUNC(tid, key).second
#define FIND_AND_CHECK_SUCCESS ds->FIND_FUNC(tid, key)
#define RQ_AND_CHECK_SUCCESS(rqcnt)                             \
  rqcnt = ds->RQ_FUNC(tid, key, key + RQSIZE - 1, rqResultKeys, \
                      (VALUE_TYPE *)rqResultValues)
#define RQ_GARBAGE(rqcnt) rqResultKeys[0] + rqResultKeys[rqcnt - 1]
#define INIT_THREAD(tid) \
  ds->initThread(tid);   \
  urcu::registerThread(tid);
#define DEINIT_THREAD(tid) \
  ds->deinitThread(tid);   \
  urcu::unregisterThread();
#define INIT_ALL urcu::init(TOTAL_THREADS);
#define DEINIT_ALL urcu::deinit(TOTAL_THREADS);

#define PRINT_OBJ_SIZES \
  cout << "sizes: node=" << (sizeof(node_t<test_type, test_type>)) << endl;

#elif defined(LAZYLIST)
#include "record_manager.h"
#include "lazylist_impl.h"

#define DS_DECLARATION lazylist<test_type, test_type, MEMMGMT_T>
#define MEMMGMT_T \
  record_manager<RECLAIM, ALLOC, POOL, node_t<test_type, test_type>>
#define DS_CONSTRUCTOR \
  new DS_DECLARATION(TOTAL_THREADS, KEY_MIN, KEY_MAX, NO_VALUE)

#define INSERT_AND_CHECK_SUCCESS \
  ds->INSERT_FUNC(tid, key, VALUE) == ds->NO_VALUE
#define DELETE_AND_CHECK_SUCCESS ds->ERASE_FUNC(tid, key) != ds->NO_VALUE
#define FIND_AND_CHECK_SUCCESS ds->FIND_FUNC(tid, key)
#define RQ_AND_CHECK_SUCCESS(rqcnt)                              \
  (rqcnt = ds->RQ_FUNC(tid, key, key + RQSIZE - 1, rqResultKeys, \
                       (VALUE_TYPE *)rqResultValues))
#define RQ_GARBAGE(rqcnt) rqResultKeys[0] + rqResultKeys[rqcnt - 1]
#define INIT_THREAD(tid) ds->initThread(tid)
#define DEINIT_THREAD(tid) ds->deinitThread(tid);
#define INIT_ALL
#define DEINIT_ALL

#define PRINT_OBJ_SIZES \
  cout << "sizes: node=" << (sizeof(node_t<test_type, test_type>)) << endl;

#elif defined(SKIPLISTLOCK)
#include "record_manager.h"
#include "skiplist_lock_impl.h"

#define DS_DECLARATION skiplist<test_type, test_type, MEMMGMT_T>
#define MEMMGMT_T                      \
  record_manager<RECLAIM, ALLOC, POOL, \
                 node_t<test_type, test_type> RQ_SNAPCOLLECTOR_OBJECT_TYPES>
#define DS_CONSTRUCTOR \
  new DS_DECLARATION(TOTAL_THREADS, KEY_MIN, KEY_MAX, NO_VALUE, glob.rngs)

#define INSERT_AND_CHECK_SUCCESS \
  ds->INSERT_FUNC(tid, key, VALUE) == ds->NO_VALUE
#define DELETE_AND_CHECK_SUCCESS ds->ERASE_FUNC(tid, key) != ds->NO_VALUE
#define FIND_AND_CHECK_SUCCESS ds->FIND_FUNC(tid, key)
#define RQ_AND_CHECK_SUCCESS(rqcnt)                              \
  (rqcnt = ds->RQ_FUNC(tid, key, key + RQSIZE - 1, rqResultKeys, \
                       (VALUE_TYPE *)rqResultValues))
#define RQ_GARBAGE(rqcnt) rqResultKeys[0] + rqResultKeys[rqcnt - 1]
#define INIT_THREAD(tid) ds->initThread(tid)
#define DEINIT_THREAD(tid) ds->deinitThread(tid);
#define INIT_ALL
#define DEINIT_ALL

#define PRINT_OBJ_SIZES                                                    \
  cout << "sizes: node="                                                   \
       << (sizeof(node_t<test_type, test_type>))RQ_SNAPCOLLECTOR_OBJ_SIZES \
       << endl;

#elif defined(LFLIST)
#include "lockfree_list_impl.h"
#include "record_manager.h"

#define DS_DECLARATION lflist<test_type, test_type, MEMMGMT_T>
#define MEMMGMT_T                      \
  record_manager<RECLAIM, ALLOC, POOL, \
                 node_t<test_type, test_type> RQ_SNAPCOLLECTOR_OBJECT_TYPES>
#define DS_CONSTRUCTOR \
  new DS_DECLARATION(TOTAL_THREADS, KEY_MIN, KEY_MAX, NO_VALUE)

#define INSERT_AND_CHECK_SUCCESS \
  ds->INSERT_FUNC(tid, key, VALUE) == ds->NO_VALUE
#define DELETE_AND_CHECK_SUCCESS ds->ERASE_FUNC(tid, key) != ds->NO_VALUE
#define FIND_AND_CHECK_SUCCESS ds->FIND_FUNC(tid, key)
#define RQ_AND_CHECK_SUCCESS(rqcnt)                             \
  rqcnt = ds->RQ_FUNC(tid, key, key + RQSIZE - 1, rqResultKeys, \
                      (VALUE_TYPE *)rqResultValues)
#define RQ_GARBAGE(rqcnt) rqResultKeys[0] + rqResultKeys[rqcnt - 1]
#define INIT_THREAD(tid) ds->initThread(tid)
#define DEINIT_THREAD(tid) ds->deinitThread(tid);
#define INIT_ALL
#define DEINIT_ALL

#define PRINT_OBJ_SIZES                                                    \
  cout << "sizes: node="                                                   \
       << (sizeof(node_t<test_type, test_type>))RQ_SNAPCOLLECTOR_OBJ_SIZES \
       << endl;

#elif defined(LFSKIPLIST)
#include "lockfree_skiplist_impl.h"
#include "record_manager.h"

#define DS_DECLARATION lfskiplist<test_type, test_type, MEMMGMT_T>
#define MEMMGMT_T                      \
  record_manager<RECLAIM, ALLOC, POOL, \
                 node_t<test_type, test_type> RQ_SNAPCOLLECTOR_OBJECT_TYPES>
#define DS_CONSTRUCTOR \
  new DS_DECLARATION(TOTAL_THREADS, KEY_MIN, KEY_MAX, NO_VALUE)

#define INSERT_AND_CHECK_SUCCESS \
  ds->INSERT_FUNC(tid, key, VALUE) == ds->NO_VALUE
#define DELETE_AND_CHECK_SUCCESS ds->ERASE_FUNC(tid, key) != ds->NO_VALUE
#define FIND_AND_CHECK_SUCCESS ds->FIND_FUNC(tid, key)
#define RQ_AND_CHECK_SUCCESS(rqcnt)                             \
  rqcnt = ds->RQ_FUNC(tid, key, key + RQSIZE - 1, rqResultKeys, \
                      (VALUE_TYPE *)rqResultValues)
#define RQ_GARBAGE(rqcnt) rqResultKeys[0] + rqResultKeys[rqcnt - 1]
#define INIT_THREAD(tid) ds->initThread(tid)
#define DEINIT_THREAD(tid) ds->deinitThread(tid);
#define INIT_ALL
#define DEINIT_ALL

#define PRINT_OBJ_SIZES                                                    \
  cout << "sizes: node="                                                   \
       << (sizeof(node_t<test_type, test_type>))RQ_SNAPCOLLECTOR_OBJ_SIZES \
       << endl;

#elif defined(RLU_LIST)
#include "record_manager.h"
#include "rlu.h"
#include "rlu_list_impl.h"

#define DS_DECLARATION rlulist<test_type, test_type>
#define MEMMGMT_T \
  record_manager<RECLAIM, ALLOC, POOL, node_t<test_type, test_type>>
#define DS_CONSTRUCTOR \
  new DS_DECLARATION(TOTAL_THREADS, KEY_MIN, KEY_MAX, NO_VALUE)

#define INSERT_AND_CHECK_SUCCESS \
  ds->INSERT_FUNC(tid, key, VALUE) == ds->NO_VALUE
#define DELETE_AND_CHECK_SUCCESS ds->ERASE_FUNC(tid, key).second
#define FIND_AND_CHECK_SUCCESS ds->FIND_FUNC(tid, key)
#define RQ_AND_CHECK_SUCCESS(rqcnt)                             \
  rqcnt = ds->RQ_FUNC(tid, key, key + RQSIZE - 1, rqResultKeys, \
                      (VALUE_TYPE *)rqResultValues)
#define RQ_GARBAGE(rqcnt) rqResultKeys[0] + rqResultKeys[rqcnt - 1]
__thread rlu_thread_data_t *rlu_self;
rlu_thread_data_t *rlu_tdata = NULL;
#define INIT_THREAD(tid)      \
  rlu_self = &rlu_tdata[tid]; \
  RLU_THREAD_INIT(rlu_self);
#define DEINIT_THREAD(tid) RLU_THREAD_FINISH(rlu_self);
#define INIT_ALL                                   \
  rlu_tdata = new rlu_thread_data_t[MAX_TID_POW2]; \
  RLU_INIT(RLU_TYPE_FINE_GRAINED, 1)
#define DEINIT_ALL \
  RLU_FINISH();    \
  delete[] rlu_tdata;

#define PRINT_OBJ_SIZES                                                  \
  cout << "sizes: node="                                                 \
       << ((sizeof(node_t<test_type, test_type>)) + RLU_OBJ_HEADER_SIZE) \
       << " including header=" << RLU_OBJ_HEADER_SIZE << endl;

#elif defined(RLU_CITRUS)
#include "record_manager.h"
#include "rlu.h"
#include "rlu_citrus_impl.h"

#define DS_DECLARATION rlucitrus<test_type, test_type>
#define MEMMGMT_T \
  record_manager<RECLAIM, ALLOC, POOL, node_t<test_type, test_type>>
#define DS_CONSTRUCTOR new DS_DECLARATION(TOTAL_THREADS, KEY_MAX, NO_VALUE)

#define INSERT_AND_CHECK_SUCCESS \
  ds->INSERT_FUNC(tid, key, VALUE) == ds->NO_VALUE
#define DELETE_AND_CHECK_SUCCESS ds->ERASE_FUNC(tid, key).second
#define FIND_AND_CHECK_SUCCESS ds->FIND_FUNC(tid, key)
#define RQ_AND_CHECK_SUCCESS(rqcnt)                             \
  rqcnt = ds->RQ_FUNC(tid, key, key + RQSIZE - 1, rqResultKeys, \
                      (VALUE_TYPE *)rqResultValues)
#define RQ_GARBAGE(rqcnt) rqResultKeys[0] + rqResultKeys[rqcnt - 1]
__thread rlu_thread_data_t *rlu_self;
rlu_thread_data_t *rlu_tdata = NULL;
#define INIT_THREAD(tid)      \
  rlu_self = &rlu_tdata[tid]; \
  RLU_THREAD_INIT(rlu_self);
#define DEINIT_THREAD(tid) RLU_THREAD_FINISH(rlu_self);
#define INIT_ALL                                   \
  rlu_tdata = new rlu_thread_data_t[MAX_TID_POW2]; \
  RLU_INIT(RLU_TYPE_FINE_GRAINED, 1)
#define DEINIT_ALL \
  RLU_FINISH();    \
  delete[] rlu_tdata;

#define PRINT_OBJ_SIZES                                                  \
  cout << "sizes: node="                                                 \
       << ((sizeof(node_t<test_type, test_type>)) + RLU_OBJ_HEADER_SIZE) \
       << " including header=" << RLU_OBJ_HEADER_SIZE << endl;

#elif defined(BUNDLE_LIST)
#define BUNDLE_TYPE_DECL LinkedBundle
#include "record_manager.h"
#include "bundle_lazylist_impl.h"

#define DS_DECLARATION bundle_lazylist<test_type, test_type, MEMMGMT_T>
#define MEMMGMT_T \
  record_manager<RECLAIM, ALLOC, POOL, node_t<test_type, test_type>>
#define DS_CONSTRUCTOR \
  new DS_DECLARATION(TOTAL_THREADS + 1, KEY_MIN, KEY_MAX, NO_VALUE)

#define INSERT_AND_CHECK_SUCCESS \
  ds->INSERT_FUNC(tid, key, VALUE) == ds->NO_VALUE
#define DELETE_AND_CHECK_SUCCESS ds->ERASE_FUNC(tid, key) != ds->NO_VALUE
#define FIND_AND_CHECK_SUCCESS ds->FIND_FUNC(tid, key)
#define RQ_AND_CHECK_SUCCESS(rqcnt)                             \
  rqcnt = ds->RQ_FUNC(tid, key, key + RQSIZE - 1, rqResultKeys, \
                      (VALUE_TYPE *)rqResultValues)
#define RQ_GARBAGE(rqcnt) rqResultKeys[0] + rqResultKeys[rqcnt - 1]
#define INIT_THREAD(tid) ds->initThread(tid)
#define INIT_RQ_THREAD(tid) ds->initThread(tid, true)
#define DEINIT_THREAD(tid) ds->deinitThread(tid);
#define DEINIT_RQ_THREAD(tid) ds->deinitThread(tid, true);
#define VALIDATE_BUNDLES                                  \
  ((DS_DECLARATION *)glob.__ds)->validateBundles(0)       \
      ? std::cout << "Bundle validation OK." << std::endl \
      : std::cout << "Bundle validation failed." << std::endl;
#define INIT_ALL
#define DEINIT_ALL VALIDATE_BUNDLES

#define BUNDLE_OBJ_SIZE (sizeof(BUNDLE_TYPE_DECL<node_t<test_type, test_type>>))
#define PRINT_OBJ_SIZES                                              \
  cout << "sizes: node="                                             \
       << ((sizeof(node_t<test_type, test_type>)) + BUNDLE_OBJ_SIZE) \
       << " including header=" << BUNDLE_OBJ_SIZE << endl;

#elif defined(BUNDLE_SKIPLIST)
#include "record_manager.h"
#include "bundle_skiplist_impl.h"

#define DS_DECLARATION bundle_skiplist<test_type, test_type, MEMMGMT_T>
#define MEMMGMT_T \
  record_manager<RECLAIM, ALLOC, POOL, node_t<test_type, test_type>>
#define DS_CONSTRUCTOR \
  new DS_DECLARATION(TOTAL_THREADS + 1, KEY_MIN, KEY_MAX, NO_VALUE, glob.rngs)

#define INSERT_AND_CHECK_SUCCESS \
  ds->INSERT_FUNC(tid, key, VALUE) == ds->NO_VALUE
#define DELETE_AND_CHECK_SUCCESS ds->ERASE_FUNC(tid, key) != ds->NO_VALUE
#define FIND_AND_CHECK_SUCCESS ds->FIND_FUNC(tid, key)
#define RQ_AND_CHECK_SUCCESS(rqcnt)                             \
  rqcnt = ds->RQ_FUNC(tid, key, key + RQSIZE - 1, rqResultKeys, \
                      (VALUE_TYPE *)rqResultValues)
#define RQ_GARBAGE(rqcnt) rqResultKeys[0] + rqResultKeys[rqcnt - 1]
#define INIT_THREAD(tid) ds->initThread(tid)
#define DEINIT_THREAD(tid) ds->deinitThread(tid);
#define VALIDATE_BUNDLES                                  \
  ((DS_DECLARATION *)glob.__ds)->validateBundles(0)       \
      ? std::cout << "Bundle validation OK." << std::endl \
      : std::cout << "Bundle validation failed." << std::endl;
#define INIT_ALL
#define DEINIT_ALL VALIDATE_BUNDLES

#define BUNDLE_OBJ_SIZE (sizeof(BUNDLE_TYPE_DECL<node_t<test_type, test_type>>))
#define PRINT_OBJ_SIZES                                              \
  cout << "sizes: node="                                             \
       << ((sizeof(node_t<test_type, test_type>)) + BUNDLE_OBJ_SIZE) \
       << " including header=" << BUNDLE_OBJ_SIZE << endl;

#elif defined(BUNDLE_CITRUS)
#define BUNDLE_TYPE_DECL LinkedBundle
#include "bundle_citrus_impl.h"
#include "record_manager.h"

#define DS_DECLARATION bundle_citrustree<test_type, test_type, MEMMGMT_T>
#define MEMMGMT_T \
  record_manager<RECLAIM, ALLOC, POOL, node_t<test_type, test_type>>
#define DS_CONSTRUCTOR new DS_DECLARATION(KEY_MAX, NO_VALUE, TOTAL_THREADS + 1)

#define INSERT_AND_CHECK_SUCCESS \
  ds->INSERT_FUNC(tid, key, VALUE) == ds->NO_VALUE
#define DELETE_AND_CHECK_SUCCESS ds->ERASE_FUNC(tid, key).second
#define FIND_AND_CHECK_SUCCESS ds->FIND_FUNC(tid, key)
#define RQ_AND_CHECK_SUCCESS(rqcnt)                             \
  rqcnt = ds->RQ_FUNC(tid, key, key + RQSIZE - 1, rqResultKeys, \
                      (VALUE_TYPE *)rqResultValues)
#define RQ_GARBAGE(rqcnt) rqResultKeys[0] + rqResultKeys[rqcnt - 1]
#define INIT_THREAD(tid) \
  ds->initThread(tid);   \
  urcu::registerThread(tid);
#define DEINIT_THREAD(tid) \
  ds->deinitThread(tid);   \
  urcu::unregisterThread();
#define VALIDATE_BUNDLES                                  \
  ((DS_DECLARATION *)glob.__ds)->validateBundles(0)       \
      ? std::cout << "Bundle validation OK." << std::endl \
      : std::cout << "Bundle validation failed." << std::endl;
#define INIT_ALL urcu::init(TOTAL_THREADS + 1);
#define DEINIT_ALL  \
  VALIDATE_BUNDLES; \
  urcu::deinit(TOTAL_THREADS + 1);

#define BUNDLE_OBJ_SIZE (sizeof(BUNDLE_TYPE_DECL<node_t<test_type, test_type>>))
#define PRINT_OBJ_SIZES                                              \
  cout << "sizes: node="                                             \
       << ((sizeof(node_t<test_type, test_type>)) + BUNDLE_OBJ_SIZE) \
       << " including header=" << BUNDLE_OBJ_SIZE << endl;

#elif defined(BUNDLE_BST)
#define BUNDLE_TYPE_DECL LinkedBundle
#define BUNDLE_LOCKFREE
#include "bundle_bst_impl.h"
#include "record_manager.h"
using namespace bundle_bst_ns;

#define DS_DECLARATION \
  bundle_bst<test_type, test_type, less<test_type>, MEMMGMT_T>
#define MEMMGMT_T \
  record_manager<RECLAIM, ALLOC, POOL, Node<test_type, test_type>>
#define DS_CONSTRUCTOR \
  new DS_DECLARATION(KEY_MAX, NO_VALUE, TOTAL_THREADS + 1, SIGQUIT)

#define INSERT_AND_CHECK_SUCCESS \
  ds->INSERT_FUNC(tid, key, VALUE) == ds->NO_VALUE
#define DELETE_AND_CHECK_SUCCESS ds->ERASE_FUNC(tid, key).second
#define FIND_AND_CHECK_SUCCESS ds->FIND_FUNC(tid, key)
#define RQ_AND_CHECK_SUCCESS(rqcnt)                               \
  (rqcnt) = ds->RQ_FUNC(tid, key, key + RQSIZE - 1, rqResultKeys, \
                        (VALUE_TYPE *)rqResultValues)
#define RQ_GARBAGE(rqcnt) rqResultKeys[0] + rqResultKeys[(rqcnt)-1]
#define INIT_THREAD(tid) ds->initThread(tid)
#define DEINIT_THREAD(tid) ds->deinitThread(tid)
#define INIT_ALL
#define DEINIT_ALL

#define BUNDLE_OBJ_SIZE (sizeof(Bundle<node_t<test_type, test_type>>))
#define PRINT_OBJ_SIZES                                          \
  cout << "sizes: node=" << (sizeof(Node<test_type, test_type>)) \
       << " descriptor=" << (sizeof(SCXRecord<test_type, test_type>)) << endl;

#elif defined(UNSAFE_LIST)
#include "unsafe_lazylist_impl.h"
#include "record_manager.h"

#define DS_DECLARATION unsafe_lazylist<test_type, test_type, MEMMGMT_T>
#define MEMMGMT_T \
  record_manager<RECLAIM, ALLOC, POOL, node_t<test_type, test_type>>
#define DS_CONSTRUCTOR \
  new DS_DECLARATION(TOTAL_THREADS, KEY_MIN, KEY_MAX, NO_VALUE)

#define INSERT_AND_CHECK_SUCCESS \
  ds->INSERT_FUNC(tid, key, VALUE) == ds->NO_VALUE
#define DELETE_AND_CHECK_SUCCESS ds->ERASE_FUNC(tid, key) != ds->NO_VALUE
#define FIND_AND_CHECK_SUCCESS ds->FIND_FUNC(tid, key)
#define RQ_AND_CHECK_SUCCESS(rqcnt)                             \
  rqcnt = ds->RQ_FUNC(tid, key, key + RQSIZE - 1, rqResultKeys, \
                      (VALUE_TYPE *)rqResultValues)
#define RQ_GARBAGE(rqcnt) rqResultKeys[0] + rqResultKeys[rqcnt - 1]
#define INIT_THREAD(tid) ds->initThread(tid)
#define INIT_RQ_THREAD(tid) ds->initThread(tid, true)
#define DEINIT_THREAD(tid) ds->deinitThread(tid);
#define DEINIT_RQ_THREAD(tid) ds->deinitThread(tid, true);
#define INIT_ALL
#define DEINIT_ALL

#define PRINT_OBJ_SIZES \
  cout << "sizes: node=" << (sizeof(node_t<test_type, test_type>)) << endl;

#elif defined(UNSAFE_SKIPLIST)
#include "record_manager.h"
#include "unsafe_skiplist_impl.h"

#define DS_DECLARATION unsafe_skiplist<test_type, test_type, MEMMGMT_T>
#define MEMMGMT_T \
  record_manager<RECLAIM, ALLOC, POOL, node_t<test_type, test_type>>
#define DS_CONSTRUCTOR \
  new DS_DECLARATION(TOTAL_THREADS, KEY_MIN, KEY_MAX, NO_VALUE, glob.rngs)

#define INSERT_AND_CHECK_SUCCESS \
  ds->INSERT_FUNC(tid, key, VALUE) == ds->NO_VALUE
#define DELETE_AND_CHECK_SUCCESS ds->ERASE_FUNC(tid, key) != ds->NO_VALUE
#define FIND_AND_CHECK_SUCCESS ds->FIND_FUNC(tid, key)
#define RQ_AND_CHECK_SUCCESS(rqcnt)                             \
  rqcnt = ds->RQ_FUNC(tid, key, key + RQSIZE - 1, rqResultKeys, \
                      (VALUE_TYPE *)rqResultValues)
#define RQ_GARBAGE(rqcnt) rqResultKeys[0] + rqResultKeys[rqcnt - 1]
#define INIT_THREAD(tid) ds->initThread(tid)
#define DEINIT_THREAD(tid) ds->deinitThread(tid);
#define INIT_ALL
#define DEINIT_ALL

#define BUNDLE_OBJ_SIZE (sizeof(Bundle<node_t<test_type, test_type>>))
#define PRINT_OBJ_SIZES \
  cout << "sizes: node=" << (sizeof(node_t<test_type, test_type>)) << endl

#elif defined(UNSAFE_CITRUS)
#include "record_manager.h"
#include "unsafe_citrus_impl.h"

#define DS_DECLARATION unsafe_citrustree<test_type, test_type, MEMMGMT_T>
#define MEMMGMT_T \
  record_manager<RECLAIM, ALLOC, POOL, node_t<test_type, test_type>>
#define DS_CONSTRUCTOR new DS_DECLARATION(KEY_MAX, NO_VALUE, TOTAL_THREADS)

#define INSERT_AND_CHECK_SUCCESS \
  ds->INSERT_FUNC(tid, key, VALUE) == ds->NO_VALUE
#define DELETE_AND_CHECK_SUCCESS ds->ERASE_FUNC(tid, key).second
#define FIND_AND_CHECK_SUCCESS ds->FIND_FUNC(tid, key)
#define RQ_AND_CHECK_SUCCESS(rqcnt)                             \
  rqcnt = ds->RQ_FUNC(tid, key, key + RQSIZE - 1, rqResultKeys, \
                      (VALUE_TYPE *)rqResultValues)
#define RQ_GARBAGE(rqcnt) rqResultKeys[0] + rqResultKeys[rqcnt - 1]
#define INIT_THREAD(tid) \
  ds->initThread(tid);   \
  urcu::registerThread(tid);
#define DEINIT_THREAD(tid) \
  ds->deinitThread(tid);   \
  urcu::unregisterThread();
#define INIT_ALL urcu::init(TOTAL_THREADS + 1);
#define DEINIT_ALL urcu::deinit(TOTAL_THREADS + 1);

#define PRINT_OBJ_SIZES \
  cout << "sizes: node=" << (sizeof(node_t<test_type, test_type>)) << endl;

#elif defined(VCASBST)
#include "record_manager.h"
#include "vcas_bst_impl.h"
using namespace vcas_bst_ns;

#define DS_DECLARATION \
  vcas_bst<test_type, test_type, less<test_type>, MEMMGMT_T>
#define MEMMGMT_T \
  record_manager<RECLAIM, ALLOC, POOL, Node<test_type, test_type>>
#define DS_CONSTRUCTOR \
  new DS_DECLARATION(KEY_MAX, NO_VALUE, TOTAL_THREADS, SIGQUIT)

#define INSERT_AND_CHECK_SUCCESS \
  ds->INSERT_FUNC(tid, key, VALUE) == ds->NO_VALUE
#define DELETE_AND_CHECK_SUCCESS ds->ERASE_FUNC(tid, key).second
#define FIND_AND_CHECK_SUCCESS ds->FIND_FUNC(tid, key)
#define RQ_AND_CHECK_SUCCESS(rqcnt)                               \
  (rqcnt) = ds->RQ_FUNC(tid, key, key + RQSIZE - 1, rqResultKeys, \
                        (VALUE_TYPE *)rqResultValues)
#define RQ_GARBAGE(rqcnt) rqResultKeys[0] + rqResultKeys[(rqcnt)-1]
#define INIT_THREAD(tid) ds->initThread(tid)
#define DEINIT_THREAD(tid) ds->deinitThread(tid)
#define INIT_ALL
#define DEINIT_ALL

#define PRINT_OBJ_SIZES                                          \
  cout << "sizes: node=" << (sizeof(Node<test_type, test_type>)) \
       << " descriptor=" << (sizeof(SCXRecord<test_type, test_type>)) << endl;

#elif defined(VCAS_LAZYLIST)

#define NVCAS_OPTIMIZATION

#include "record_manager.h"
#include "vcas_lazylist_impl.h"
using namespace vcas_lazylist;

#define DS_DECLARATION \
  lazylist<test_type, test_type, MEMMGMT_T>
#define MEMMGMT_T \
  record_manager<RECLAIM, ALLOC, POOL, node_t<test_type, test_type> >
#define DS_CONSTRUCTOR \
  new DS_DECLARATION(TOTAL_THREADS, KEY_MIN, KEY_MAX, NO_VALUE)

#define INSERT_AND_CHECK_SUCCESS \
  ds->INSERT_FUNC(tid, key, VALUE) == ds->NO_VALUE
#define DELETE_AND_CHECK_SUCCESS ds->ERASE_FUNC(tid, key) != ds->NO_VALUE
#define FIND_AND_CHECK_SUCCESS ds->FIND_FUNC(tid, key)
#define RQ_AND_CHECK_SUCCESS(rqcnt)                              \
  (rqcnt = ds->RQ_FUNC(tid, key, key + RQSIZE - 1, rqResultKeys, \
                       (VALUE_TYPE *)rqResultValues))
#define RQ_GARBAGE(rqcnt) rqResultKeys[0] + rqResultKeys[(rqcnt)-1]
#define INIT_THREAD(tid) ds->initThread(tid)
#define DEINIT_THREAD(tid) ds->deinitThread(tid)
#define INIT_ALL
#define DEINIT_ALL

#define PRINT_OBJ_SIZES \
  cout << "sizes: node=" << (sizeof(node_t<test_type, test_type>)) << endl;

#elif defined(VCAS_SKIPLIST)

#define NVCAS_OPTIMIZATION

#include "record_manager.h"
#include "vcas_skiplist_lock_impl.h"

using namespace vcas_skiplist_lock;

#define DS_DECLARATION skiplist<test_type, test_type, MEMMGMT_T>
#define MEMMGMT_T record_manager<RECLAIM, ALLOC, POOL, node_t<test_type, test_type> RQ_SNAPCOLLECTOR_OBJECT_TYPES>
#define DS_CONSTRUCTOR new DS_DECLARATION(TOTAL_THREADS, KEY_MIN, KEY_MAX, NO_VALUE, glob.rngs)
#define INSERT_AND_CHECK_SUCCESS ds->INSERT_FUNC(tid, key, VALUE) == ds->NO_VALUE
#define DELETE_AND_CHECK_SUCCESS ds->ERASE_FUNC(tid, key) != ds->NO_VALUE
#define FIND_AND_CHECK_SUCCESS ds->FIND_FUNC(tid, key)
#define RQ_AND_CHECK_SUCCESS(rqcnt)                              \
  (rqcnt = ds->RQ_FUNC(tid, key, key + RQSIZE - 1, rqResultKeys, \
                       (VALUE_TYPE *)rqResultValues))
#define RQ_GARBAGE(rqcnt) rqResultKeys[0] + rqResultKeys[rqcnt - 1]
#define INIT_THREAD(tid) ds->initThread(tid)
#define DEINIT_THREAD(tid) ds->deinitThread(tid);
#define INIT_ALL
#define DEINIT_ALL

#define PRINT_OBJ_SIZES                                                    \
  cout << "sizes: node="                                                   \
       << (sizeof(node_t<test_type, test_type>))RQ_SNAPCOLLECTOR_OBJ_SIZES \
       << endl;

///////////
///////////
/////////// TODO: what needs to change here ???
/////////// 1. no longer find, insert, contains... rather addEdge, addVertex, removeEdge, etc.
/////////// I think I just need to make a new way of filling the DS in main and new way of performing the trial
/////////// in terms of which operations are being performed
/////////// OR: directly call the methods from main ? this could be easier

#elif defined(VCAS_GRAPH)

#define NVCAS_OPTIMIZATION

#include "record_manager.h"
#include "vcas_graph.h"

using namespace vcas_graph;

#define DS_DECLARATION graph<test_type, test_type, MEMMGMT_T> // #TODO: check how the others define these things, make sure mine matches the def here
#define MEMMGMT_T record_manager<RECLAIM, ALLOC, POOL, node_t<test_type, test_type> RQ_SNAPCOLLECTOR_OBJECT_TYPES>
#define DS_CONSTRUCTOR new DS_DECLARATION(TOTAL_THREADS, KEY_MIN, KEY_MAX, NO_VALUE, glob.rngs)
#define INSERT_VERTEX_AND_CHECK_SUCCESS ds->INSERT_VERTEX_FUNC(tid, key, VALUE) == ds->NO_VALUE // TODO: make sure the functions match this signature
#define DELETE_VERTEX_AND_CHECK_SUCCESS ds->REMOVE_VERTEX_FUNC(tid, key) != ds->NO_VALUE
#define FIND_VERTEX_AND_CHECK_SUCCESS ds->FIND_VERTEX_FUNC(tid, key)

#define INSERT_EDGE_AND_CHECK_SUCCESS ds->INSERT_EDGE_FUNC(tid, key_s, key_e) == ds->NO_VALUE // TODO: make sure the functions match this signature
#define DELETE_EDGE_AND_CHECK_SUCCESS ds->REMOVE_EDGE_FUNC(tid, key_s, key_e) != ds->NO_VALUE
#define FIND_EDGE_AND_CHECK_SUCCESS ds->FIND_EDGE_FUNC(tid, key_s, key_e)

// TODO: [future todo] figure out RQs in the scope of graphs
#define RQ_AND_CHECK_SUCCESS(rqcnt) (rqcnt = ds->RQ_FUNC(tid, key, key + RQSIZE - 1, rqResultKeys, (VALUE_TYPE *)rqResultValues))
#define RQ_GARBAGE(rqcnt) rqResultKeys[0] + rqResultKeys[rqcnt - 1]
#define INIT_THREAD(tid) ds->initThread(tid)
#define DEINIT_THREAD(tid) ds->deinitThread(tid);
#define INIT_ALL
#define DEINIT_ALL

#define PRINT_OBJ_SIZES                                                    \
  cout << "sizes: node="                                                   \
       << (sizeof(node_t<test_type, test_type>))RQ_SNAPCOLLECTOR_OBJ_SIZES \
       << endl;

///////////
///////////
///////////

#elif defined(VCAS_CITRUS)

#define NVCAS_OPTIMIZATION

#include "record_manager.h"
#include "vcas_citrus_impl.h"

using namespace vcas_citrus;

#define DS_DECLARATION citrustree<test_type, test_type, MEMMGMT_T>
#define MEMMGMT_T \
  record_manager<RECLAIM, ALLOC, POOL, node_t<test_type, test_type>>
#define DS_CONSTRUCTOR new DS_DECLARATION(MAXKEY, NO_VALUE, TOTAL_THREADS)

#define INSERT_AND_CHECK_SUCCESS \
  ds->INSERT_FUNC(tid, key, VALUE) == ds->NO_VALUE
#define DELETE_AND_CHECK_SUCCESS ds->ERASE_FUNC(tid, key).second
#define FIND_AND_CHECK_SUCCESS ds->FIND_FUNC(tid, key)
#define RQ_AND_CHECK_SUCCESS(rqcnt)                             \
  rqcnt = ds->RQ_FUNC(tid, key, key + RQSIZE - 1, rqResultKeys, \
                      (VALUE_TYPE *)rqResultValues)
#define RQ_GARBAGE(rqcnt) rqResultKeys[0] + rqResultKeys[rqcnt - 1]
#define INIT_THREAD(tid) \
  ds->initThread(tid);   \
  urcu::registerThread(tid);
#define DEINIT_THREAD(tid) \
  ds->deinitThread(tid);   \
  urcu::unregisterThread();
#define INIT_ALL urcu::init(TOTAL_THREADS);
#define DEINIT_ALL urcu::deinit(TOTAL_THREADS);

#define PRINT_OBJ_SIZES \
  cout << "sizes: node=" << (sizeof(node_t<test_type, test_type>)) << endl;
#else
#error "Failed to define a data structure"
#endif

#endif /* DATA_STRUCTURE_H */
