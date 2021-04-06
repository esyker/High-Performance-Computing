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
    int* distances,
    bool* visited)
{

    #pragma omp parallel for schedule(dynamic, 64)
    for (int i=0; i<frontier->count; i++) {

        int node = frontier->vertices[i];

        int start_edge = g->outgoing_starts[node];
        int end_edge = (node == g->num_nodes - 1)
                           ? g->num_edges
                           : g->outgoing_starts[node + 1];

        int new_frontier_prev[end_edge-start_edge];
        int count = 0;

        // attempt to add all neighbors to the new frontier
        for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
            int outgoing = g->outgoing_edges[neighbor];

            if(visited[outgoing] || distances[outgoing] != NOT_VISITED_MARKER)continue;
            if(__sync_bool_compare_and_swap(&distances[outgoing], NOT_VISITED_MARKER, distances[node] + 1)){
                new_frontier_prev[count] = outgoing;
                count += 1;
            }
        }
        if(count > 0) {
            int index = __sync_fetch_and_add(&new_frontier->count, count);
            for (int j=0; j<count; j++)
                new_frontier->vertices[index+j] = new_frontier_prev[j];
        }
    }

    #pragma omp parallel for
    for(int i = 0; i < new_frontier->count; i++)
        visited[new_frontier->vertices[i]] = true;
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
    bool* visited = (bool*) malloc(graph->num_nodes * sizeof(bool));
    if(visited == NULL){
        printf("Error allocating memory for visited array!\n");
        exit(-1);
    }


    // initialize all nodes to NOT_VISITED
    #pragma omp parallel for
    for (int i=0; i<graph->num_nodes; i++){
        sol->distances[i] = NOT_VISITED_MARKER;
        visited[i] = false;
    }

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);

        if(frontier->count >= 64)
            top_down_step_parallel(graph, frontier, new_frontier, sol->distances, visited);
        else
            top_down_step(graph, frontier, new_frontier, sol->distances);

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

/*
for each vertex v in graph:
    if v has not been visited AND 
       v shares an incoming edge with a vertex u on the frontier:
          add vertex v to frontier;
*/
void bottom_up_step(
    graph* g,
    vertex_set* frontier,
	vertex_set* new_frontier,
    int* distances) {
    
	for (int i=0; i < g->num_nodes; i++) {                   
		if (distances[i] == NOT_VISITED_MARKER) {
			int start_edge = g->incoming_starts[i];
			int end_edge = (i == g->num_nodes-1)? g->num_edges : g->incoming_starts[i + 1];
			for(int neighbor = start_edge; neighbor < end_edge; neighbor++) {
				for (int j=0;j<frontier->count;j++){
					if(frontier->vertices[j]==neighbor){
						distances[i]=distances[neighbor]-1;
						int index=new_frontier->count++;
						new_frontier->vertices[index]=i;
					}
				}
				/*
				int incoming = g->incoming_edges[neighbor];
				//int node = frontier->present[i];
				//if(__sync_bool_compare_and_swap(&frontier->present[incoming], iteration, distances[node] + 1)) {
				if(frontier->vertices[incoming] == ) {
					distances[i] = distances[incoming] + 1;                        
					//local_count ++;
					int index = new_frontier->count++;
					new_frontier->vertices[index] = ;
					break;
				}*/
			}
		}
	}
}

void bfs_bottom_up(Graph graph, solution* sol)
{
	vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);
	
    vertex_set* frontier = &list1; 
    vertex_set* new_frontier = &list2;
	
    bool* visited = (bool*) malloc(graph->num_nodes * sizeof(bool));
    if(visited == NULL){
        printf("Error allocating memory for visited array!\n");
        exit(-1);
    }
    // initialize all nodes to NOT_VISITED
    #pragma omp parallel for
    for (int i=0; i<graph->num_nodes; i++){
        sol->distances[i] = NOT_VISITED_MARKER;
        visited[i] = false;
    }
	
	frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    
    // printf("!!!!!!!!!!!!!!!!!!!!fro2: %-10d\n", frontier->count);
    // just like pop the queue
    while (frontier->count != 0) {
        vertex_set_clear(new_frontier);
        bottom_up_step(graph, frontier, new_frontier, sol->distances);
		// swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;

    }
}

void bfs_hybrid(Graph graph, solution* sol)
{
    // TODO STUDENTS:
    //
    // You will need to implement the "hybrid" BFS here as
    // described in the handout.
}
