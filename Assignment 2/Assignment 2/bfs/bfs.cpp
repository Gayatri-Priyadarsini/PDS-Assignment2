#include "bfs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstddef>
#include <omp.h>

#include "../common/CycleTimer.h"
#include "../common/graph.h"

#define ROOT_NODE_ID 0
#define NOT_VISITED_MARKER -1

void vertex_set_clear(vertex_set* list) {
    list->count = 0;
}

void vertex_set_init(vertex_set* list, int count) {
    list->max_vertices = count;
    list->vertices = (int*)malloc(sizeof(int) * list->max_vertices);
    vertex_set_clear(list);
}

// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void top_down_step(
    Graph g,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{

    for (int i=0; i<frontier->count; i++) {

        int node = frontier->vertices[i];

        int start_edge = g->outgoing_starts[node];
        int end_edge = (node == g->num_nodes - 1)
                           ? g->num_edges
                           : g->outgoing_starts[node + 1];

        // attempt to add all neighbors to the new frontier
        for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
            int outgoing = g->outgoing_edges[neighbor];

            if (distances[outgoing] == NOT_VISITED_MARKER) {
                distances[outgoing] = distances[node] + 1;
                int index = new_frontier->count++;
                new_frontier->vertices[index] = outgoing;
            }
        }
    }
}

void top_down_step_parallel(
    Graph g,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{   int new_frontier_distance = distances[frontier->vertices[0]] + 1;
    
    #pragma omp parallel 
    {

        int tid=omp_get_thread_num();
        int nthreads=omp_get_num_threads();
        for (int i=tid; i<frontier->count; i+=nthreads) {

            int node = frontier->vertices[i];
            int curr_level=distances[node];

            int start_edge = g->outgoing_starts[node];
            int end_edge = (node == g->num_nodes - 1)
                               ? g->num_edges
                               : g->outgoing_starts[node + 1];

            //idk
            int num_neighbours=end_edge-start_edge;
            int to_be_added[num_neighbours];
            int counter =0;

            // vertex_set list3;
            // vertex_set_init(&list3,num_neighbours);
            // vertex_set* local_new_frontier=&list3;

            // attempt to add all neighbors to the new frontier
            for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
                int outgoing = g->outgoing_edges[neighbor];

                if (distances[outgoing] == NOT_VISITED_MARKER) {
                    int idx;
                    if (distances[outgoing]==NOT_VISITED_MARKER)
                    { 
                        #pragma omp atomic capture
                        {
                            idx=distances[outgoing];
                            distances[outgoing]=new_frontier_distance;
                        }
                    }
                    if (idx==NOT_VISITED_MARKER)
                    {
                        to_be_added[counter] = outgoing;
                        counter++;
                    }
                }
            }

            // #pragma omp critical
            // {   
            //     for (int j=0; j<local_new_frontier->count; j++) {
            //         new_frontier->vertices[new_frontier->count+j] = local_new_frontier->vertices[j];
            //         new_frontier->count++;
            //     }
            // }
            if (counter > 0) {
                //int index = __sync_fetch_and_add(&new_frontier->count, counter);
                int index;
                #pragma omp atomic capture
                {   
                    index=new_frontier->count;
                    new_frontier->count+=counter;    
                }
            for (int j=0; j<counter; j++) {
                new_frontier->vertices[index+j] = to_be_added[j];
            }
            }
        }
    }
}

// Implements top-down BFS.
//
// Result of execution is that, for each node in the graph, the
// distance to the root is stored in sol.distances.
void bfs_top_down(Graph graph, solution* sol) {

    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);

        top_down_step_parallel(graph, frontier, new_frontier, sol->distances);

#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }
    //free(visited);
}

void update_visited(vertex_set* frontier, bool* visited) {
    for (size_t i = 0; i < frontier->count; i++) {
        int index = frontier->vertices[i];
        visited[index] = true;
    }
}

void bottom_up_step(
    Graph g,
    bool* visited,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{
    // obtain the distance for the current frontier
    int new_frontier_distance = distances[frontier->vertices[0]] + 1;

    // loop through all the vertices in the graph
    for (int i = 0; i < g->num_nodes; i++) {
        if (visited[i]) continue;
        // if (distances[i]==0) continue;
        int start_edge = g->incoming_starts[i];
        int end_edge = (i == g->num_nodes - 1)
                        ? g->num_edges
                        : g->incoming_starts[i + 1];
        // loop through all the neighbours for the given vertex
        for (int neighbour = start_edge; neighbour < end_edge; neighbour++) {
            int incoming = g->incoming_edges[neighbour];
            // check if any neighbour has already been visited
            // if true, this means that the neighbour is on the frontier
            // and i should be added to the new_frontier
            if (!visited[incoming]) continue;
            // if (distances[incoming]==NOT_VISITED_MARKER) continue;
            int index = new_frontier->count++;
            new_frontier->vertices[index] = i;
            distances[i] = new_frontier_distance;
            break;
        }
    }
}

void bottom_up_step_parallel(
    Graph g,
    bool* visited,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{

    // obtain the distance for the current frontier
    int new_frontier_distance = distances[frontier->vertices[0]] + 1;

    int nodes_per_thread = 256; 

    #pragma omp parallel for schedule(dynamic, 1)
    for (int j = 0; j < g->num_nodes; j += nodes_per_thread) {
        // for the last iteration, if the number of nodes assigned
        // exceeds the maximum, set a ceiling to the size
        int size = (nodes_per_thread + j >= g->num_nodes) ? g->num_nodes - j : nodes_per_thread;
        int counter = 0;
        int to_be_added[size];

        for (int i = j; i < j + size; i++) {
            if (visited[i]) continue;
            int start_edge = g->incoming_starts[i];
            int end_edge = (i == g->num_nodes - 1)
                ? g->num_edges
                : g->incoming_starts[i + 1];
            for (int neighbour=start_edge; neighbour<end_edge; neighbour++) {
                int incoming = g->incoming_edges[neighbour];
                // check if any neighbour has already been visited
                // if true, this means that the neighbour is on the frontier
                // and i should be added to the new_frontier
                if (!visited[incoming]) continue;
                to_be_added[counter] = i;
                counter++;
                break;
            }
        }

        //int index = __sync_fetch_and_add(&new_frontier->count, counter);
        int index;
        #pragma omp atomic capture
        {   
            index=new_frontier->count;
            new_frontier->count+=counter;    
        }
        for (int i = 0; i < counter; i++){
            new_frontier->vertices[index + i] = to_be_added[i];
            distances[to_be_added[i]] = new_frontier_distance;
        }
    }
}



void bfs_bottom_up(Graph graph, solution* sol)
{
    // You will need to implement the "bottom up" BFS here as
    // described in the handout.
    //
    // As a result of your code's execution, sol.distances should be
    // correctly populated for all nodes in the graph.
    //
    // As was done in the top-down case, you may wish to organize your
    // code by creating subroutine bottom_up_step() that is called in
    // each step of the BFS process.

    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;

    // create a new array that shows if a node has been visited
    // which is only updated at the end of each iteration
    bool* visited = (bool*) malloc(graph->num_nodes * sizeof(bool));

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++)
    {
        sol->distances[i] = NOT_VISITED_MARKER;
        visited[i] = false;
    }

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    visited[ROOT_NODE_ID] = true;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);

        bottom_up_step_parallel(graph, visited,frontier, new_frontier, sol->distances);

        update_visited(new_frontier, visited);

#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }
    free(visited);
}

void bfs_hybrid(Graph graph, solution* sol)
{
    // You will need to implement the "hybrid" BFS here as
    // described in the handout.

    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;

    // create a new array that shows if a node has been visited
    // which is only updated at the end of each iteration
    bool* visited = (bool*) malloc(graph->num_nodes * sizeof(bool));

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++)
    {
        sol->distances[i] = NOT_VISITED_MARKER;
        visited[i] = false;
    }

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    visited[ROOT_NODE_ID] = true;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0) {

        //------------------------------------

        //     //idk
        // int m_f=0;
        // int m_u=0;
        // int n_f=frontier->count;
        // int alpha=14;
        // int beta=24;

        // for(int i=0;i<frontier->count;i++)
        // {
        //     int start_edge = graph->outgoing_starts[i];
        //     int end_edge = (i == graph->num_nodes - 1)
        //                        ? graph->num_edges
        //                        : graph->outgoing_starts[i + 1];

        //     int num_neighbours=end_edge-start_edge;
        //     m_f=m_f+num_neighbours;
        // }

        // for(int i=0;i<graph->num_nodes;i++)
        // {   
        //     if(sol->distances[i]==NOT_VISITED_MARKER)
        //     {
        //         int start_edge = graph->outgoing_starts[i];
        //         int end_edge = (i == graph->num_nodes - 1)
        //                            ? graph->num_edges
        //                            : graph->outgoing_starts[i + 1];

        //         int num_neighbours=end_edge-start_edge;
        //         m_u=m_u+num_neighbours;
        //     }
        // }

        //printf("%d,%d",m_f,m_u);

        // int c_tb= m_u/alpha;
        // int c_bt= (graph->num_nodes)/beta;

        // printf("%d,%d",c_tb,c_bt);

        //----------------------------------------

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);

        // if (m_f>c_tb) {
        //     //printf("bottom_up");
        //     bottom_up_step_parallel(graph, visited, frontier, new_frontier, sol->distances);
        // } else if(n_f <c_bt){
        //     //printf("top_down");
        //     top_down_step_parallel(graph, frontier, new_frontier, sol->distances);
        // }

        if ((double)frontier->count / graph->num_nodes > 0.03) {
            bottom_up_step_parallel(graph, visited, frontier, new_frontier, sol->distances);
        } else{
            top_down_step_parallel(graph, frontier, new_frontier, sol->distances);
        }
        update_visited(new_frontier, visited);

#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }
    free(visited);
}
