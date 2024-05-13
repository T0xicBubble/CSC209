#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "streets.h"
#define MAX 999999

struct node {
    int id;
    int num_ways;
    double lat;
    double lon;
    int * way_ids;
};

struct way {
    int id;
    int num_nodes;
    char * name;
    float maxspeed;
    bool oneway;
    int * node_ids;
};

struct ssmap {
    int num_nodes;
    int num_ways;
    struct node ** nodes;
    struct way ** ways;
    // Used to insert nodes and ways into their corresponding arrays
    int curr_node_pos;
    int curr_way_pos;
};

// Used for priority queue implementation and create_path
struct vertex {
    double dist;
    bool visited;
    struct node *node;
    struct vertex *prev;
};


// Creates and returns a ssmap, param information in streets.h
struct ssmap * 
ssmap_create(int nr_nodes, int nr_ways)
{
    struct ssmap * new_map = malloc(sizeof(struct ssmap));

    // Check for malloc error
    if (new_map == NULL) 
    {
        free(new_map);
        return NULL;
    }

    // Intialize map
    new_map->num_ways = nr_ways;
    new_map->num_nodes = nr_nodes;
    new_map->curr_way_pos = 0;
    new_map->curr_node_pos = 0;
    new_map->nodes = NULL;
    new_map->ways = NULL;

    // Initialize space for the nodes and ways
    if (nr_nodes > 0 && nr_ways > 0) 
    {
        new_map->nodes = malloc(sizeof(struct node *) * nr_nodes);
        new_map->ways = malloc(sizeof(struct way *) * nr_ways);
    }

    // Check for malloc error
    if (new_map->nodes == NULL || new_map->ways == NULL) 
    {
        free(new_map->nodes);
        free(new_map->ways);
        return NULL;
    }

    return new_map;
}


// Not needed
bool
ssmap_initialize(struct ssmap * m)
{ 
    return true;
}


// Frees all nodes and ways in map and then the map itself
void
ssmap_destroy(struct ssmap * m)
{
    // Free nodes
    for (int i = 0; i < m->num_nodes; i++) 
    {
        if(m->nodes[i] != NULL) {
            free(m->nodes[i]->way_ids);
            free(m->nodes[i]);
        }
    }

    // Free ways
    for (int i = 0; i < m->num_ways; i++) 
    {
        if (m->ways[i] != NULL) 
        {
            free(m->ways[i]->node_ids);
            free(m->ways[i]->name);
            free(m->ways[i]);
        }
    }
       
    // Free map 
    free(m->ways);
    free(m->nodes);
    free(m);
}


// Creates and adds a way to map and returns said way
struct way * 
ssmap_add_way(struct ssmap * m, int id, const char * name, float maxspeed, bool oneway, 
              int num_nodes, const int node_ids[num_nodes])
{
    struct way * new_way = malloc(sizeof(struct way));

    // Check for malloc error
    if (new_way == NULL) 
    {
        free(new_way);
        return NULL;
    }

    new_way->id = id;
    new_way->num_nodes = num_nodes;
    new_way->name = malloc(sizeof(char) * strlen(name) + 1);
    
    // Check for malloc error
    if (new_way->name == NULL) {

        free(new_way->name);
        free(new_way);
        return NULL;
    }

    strncpy(new_way->name, name, sizeof(char) * strlen(name) + 1);
    new_way->maxspeed = maxspeed;
    new_way->oneway = oneway;
    new_way->node_ids = malloc(sizeof(int) * num_nodes);
    
    // Check for malloc error
    if (new_way->node_ids == NULL) 
    {
        free(new_way->node_ids);
        free(new_way);
        return NULL;
    }

    // Add associated nodes
    for (int i = 0; i < num_nodes; i++) 
    {
        new_way->node_ids[i] = node_ids[i];
    }

    // Add way to map
    m->ways[m->curr_way_pos] = new_way;
    m->curr_way_pos++;
    return new_way;
}


// Creates and adds a node to map and returns said node
struct node * 
ssmap_add_node(struct ssmap * m, int id, double lat, double lon, 
               int num_ways, const int way_ids[num_ways])
{
    struct node * new_node = malloc(sizeof(struct node));
    
    // Check for malloc error
    if (new_node == NULL) 
    {
        free(new_node);
        return NULL;
    }

    new_node->id = id;
    new_node->num_ways = num_ways;
    new_node->lat = lat;
    new_node->lon = lon;
    new_node->way_ids = malloc(sizeof(int) * num_ways);
    
    // Check for malloc error
    if (new_node->way_ids == NULL) 
    {
        free(new_node->way_ids);
        free(new_node);
        return NULL;
    }

    // Add associated ways
    for (int i = 0; i < num_ways; i++) 
    {
        new_node->way_ids[i] = way_ids[i];
    }

    // Add node to map
    m->nodes[m->curr_node_pos] = new_node;
    m->curr_node_pos += 1;
    return new_node;
}


/**
 * Finds the way corresponding to a given way_id in map
 *
 * @param m The ssmap structure where the way object should be located.
 * @param way_id The id of the way object to be found.
 * @return The way corresponding to the given way_id, if the way is not found returns NULL.
 */
struct way * 
find_way(const struct ssmap * m, const int way_id) 
{
    struct way * found_way = NULL;
    for (int i = 0; i < m->num_ways; i++) 
    {
        if (m->ways[i]->id == way_id) 
        {
            found_way = m->ways[i];
        }
    }
    if (found_way == NULL) {
        return NULL;
    }

    return found_way;
}


/**
 * Finds the node corresponding to a given node_id in map
 *
 * @param m The ssmap structure where the way object should be located.
 * @param node_id The id of the way object to be found.
 * @return The way corresponding to the given node_id, if the node is not found returns NULL.
 */
struct node * 
find_node(const struct ssmap * m, const int node_id) 
{
    struct node * found_node = NULL;
    for (int i = 0; i < m->num_nodes; i++) 
    {
        if (m->nodes[i]->id == node_id) 
        {
            found_node = m->nodes[i];
        }
    }
    if (found_node == NULL) 
    {
        return NULL;
    }

    return found_node;
}


// Prints the way corresponding to given <id>
void
ssmap_print_way(const struct ssmap * m, int id)
{
    struct way * found_way = find_way(m, id);

    if (found_way != NULL) 
    {
        printf("Way %d: %s\n", id, found_way->name);
    } 
    else
    {
        printf("error: way %d does not exist.\n", id);
    }
}


// Prints the node corresponding to given <id>
void
ssmap_print_node(const struct ssmap * m, int id)
{
    struct node * found_node = find_node(m, id);

    if (found_node != NULL) 
    {
        printf("Node %d: (%.7lf, %.7lf)\n", id, found_node->lat, found_node->lon);
    } 
    else
    {
        printf("error: node %d does not exist.\n", id);
    }
}


// Prints all ways that contain <name> in their names
void 
ssmap_find_way_by_name(const struct ssmap * m, const char * name)
{   
    bool first = true;

    // Iterate through ways to check for instances of <name>
    for (int i = 0; i < m->num_ways; i++) 
    {
        // Use strstr to see if current way contains <name>
        if (strstr(m->ways[i]->name, name) != NULL) 
        {
            if (first) 
            {
                printf("%d", m->ways[i]->id);
                first = false;
            } 
            else 
            {
                printf(" %d", m->ways[i]->id);
            }
        } 
    }

    printf("\n");
}


/**
 * Prints out all elements of an integer array seperated by spaces.
 * 
 * @param arr The array to be printed
 * @param size The size of <arr>
 */
void
print_int_array(const int * arr, const int size) {
    bool first = true;
    
    for (int i = 0; i < size; i++) 
    {
        if (first) 
        {
            printf("%d", arr[i]);
            first = false;
        } 
        else 
        {
            printf(" %d", arr[i]);
        }
    }

    printf("\n");
}


/**
 * Checks if a given integer is in an array of integers.
 * 
 * @param haystack The array of integers to be searched in.
 * @param size The size of the haystack
 * @param needle The integer to be found
 * @return True if the needle is found and false if it is not found within the haystack
 */
bool 
find_in_int_array(const int * haystack, const int size, const int needle) 
{
    for (int i = 0; i < size; i++) 
    {
        if (haystack[i] == needle) 
        {
            return true;
        }
    }

    return false;
}


/**
 * Removes duplicate way ids given two arrays of way ids.
 * Removes the duplicates by removing the duplicate from one array only alternating, this way
 * the duplicate values still exist but are only in one array. 
 * 
 * @param arr1 The first array.
 * @param arr2 The second array
 * @param size1 The size of the first array
 * @param size2 The size of the second array
 */
void 
remove_duplicates_ways(int * arr1, int * arr2, const int size1, const int size2)
{
    bool alternate = true;

    for (int i = 0; i < size1; i++) 
    {
        for (int j = 0; j < size2; j++) 
        {
            if (arr1[i] == arr2[j]) 
            {
                if (alternate) 
                {
                    arr1[i] = -1;
                    alternate = false;
                } 
                else 
                {
                    arr2[j] = -1;
                    alternate = true;
                }
            }
        }
    }
}


/**
 * Finds all duplicate values from two integer arrays and returns the
 * duplicated values in the form of an integer array.
 * 
 * The size of the array containing duplicates is then stored in the <size1> parameter
 * 
 * @param arr1 The first array.
 * @param arr2 The second array
 * @param size1 The pointer to the size of the first array, the value which will hold the size
 * of the new array after the function returns
 * @param size2 The size of the second array
 * @return Returns an array containing all the values at appear in arr1 and arr2
 */
int * 
find_duplicates(int * arr1, int * arr2, int * size1, int size2)
{
    int max_size = 0;

    // Find the larger value between <size1> and <size2>
    if (* size1 > size2) 
    {
        max_size = * size1;
    }

    max_size = size2;

    int * duplicates = malloc(sizeof(int) * max_size);
    
    // Check for malloc error
    if (duplicates == NULL) 
    {
        free(duplicates);
        return NULL;
    }

    int count = 0;

    // Search for duplicates.
    for (int i = 0; i < * size1; i++) 
    {
        for (int j = 0; j < size2; j++) 
        {
            if (arr1[i] == arr2[j]) 
            {
                duplicates[count] = arr1[i];
                count++;
            }
        }
    }

    // Store the size of the duplicates array in <size1>
    * size1 = count;
    return duplicates;
}


/**
 * Finds all nodes associated to a given array of ways and stores it in the
 * parameter <found_nodes>
 * 
 * @param m The ssmap structure where the ways are located
 * @param ways The array of ways
 * @param way_count The size of ways
 * @param found_nodes The array to store all the associated nodes
 * @param found_node_count The value in which the size of <found_nodes> will be stored in
 */
void 
find_associated_nodes(const struct ssmap * m, int * ways, int way_count, int * found_nodes, int * found_node_count)
{
    int node_count = * found_node_count;
    struct way * curr_way;

    // Iterate through ways
    for (int i = 0; i < m->num_ways; i++) 
    {
        curr_way = m->ways[i];

        // Make sure curr_way exists in ways
        if (find_in_int_array(ways, way_count, curr_way->id)) 
        {
            // Iterate through all nodes corresponding to curr_way
            for (int j = 0; j < curr_way->num_nodes; j++) 
            {
                // Check if node already exists inside found_nodes
                if (!find_in_int_array(found_nodes, node_count, curr_way->node_ids[j])) 
                {
                    found_nodes[node_count] = curr_way->node_ids[j];
                    node_count++;
                }
            }
        } 
    }

    * found_node_count = node_count;   
}


// Finds all nodes that contain at least two ways, at least one corresponding to ways that contain <name1> 
// and at least one from a way containing <name2>.
// If <name2> is null, finds all nodes corresponding to ways that contain <name1>
void 
ssmap_find_node_by_names(const struct ssmap * m, const char * name1, const char * name2)
{
    int found_nodes[m->num_nodes];
    int node_count = 0;

    // Initialize found_nodes 
    for (int i = 0; i < m->num_nodes; i++) 
    {
        found_nodes[i] = -1;
    }

    int ways[m->num_ways];
    int way_count = 0;

    // Find assosiated ways for <name1>
    for (int i = 0; i < m->num_ways; i++) 
    {
        if (strstr(m->ways[i]->name, name1) != NULL) 
        {
            ways[way_count] = m->ways[i]->id;
            way_count++;
        }
    }

    if (name2 == NULL) 
    {
        // Find associated nodes
        find_associated_nodes(m, ways, way_count, found_nodes, &node_count);

        // Print associated nodes
        print_int_array(found_nodes, node_count);
    } 
    else 
    {
        int found_nodes2[m->num_nodes];
        int node_count2 = 0;

        // Initialize found_nodes2 
        for (int i = 0; i < m->num_nodes; i++) 
        {
            found_nodes2[i] = -1;
        }

        int ways2[m->num_ways];
        int way_count2 = 0;

        // Find assosiated ways for <name2>
        for (int i = 0; i < m->num_ways; i++) 
        {
            if (strstr(m->ways[i]->name, name2) != NULL) 
            {
                ways2[way_count2] = m->ways[i]->id;
                way_count2++;
            }
        }

        // Remove duplicate ways corresponding to <name1> and <name2>
        remove_duplicates_ways(ways, ways2, way_count, way_count2); 

        // Find associated nodes given the new ways with "removed" duplicates
        find_associated_nodes(m, ways, way_count, found_nodes, &node_count);
        find_associated_nodes(m, ways2, way_count2, found_nodes2, &node_count2);

        // Find nodes that are found in both found_nodes and found_nodes2
        int * in_both_lsts = find_duplicates(found_nodes, found_nodes2, &node_count, node_count2);
        print_int_array(in_both_lsts, node_count);
        free(in_both_lsts);
    }  
}


/**
 * Converts from degree to radian
 *
 * @param deg The angle in degrees.
 * @return the equivalent value in radian
 */
#define d2r(deg) ((deg) * M_PI/180.)


/**
 * Calculates the distance between two nodes using the Haversine formula.
 *
 * @param x The first node.
 * @param y the second node.
 * @return the distance between two nodes, in kilometre.
 */
static double
distance_between_nodes(const struct node * x, const struct node * y) 
{
    double R = 6371.;       
    double lat1 = x->lat;
    double lon1 = x->lon;
    double lat2 = y->lat;
    double lon2 = y->lon;
    double dlat = d2r(lat2-lat1); 
    double dlon = d2r(lon2-lon1); 
    double a = pow(sin(dlat/2), 2) + cos(d2r(lat1)) * cos(d2r(lat2)) * pow(sin(dlon/2), 2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a)); 
    return R * c; 
}


// Given an ordered array of node ids, returns the total travel time of said array
double 
ssmap_path_travel_time(const struct ssmap * m, int size, int node_ids[size])
{
    double total_travel_time = 0;
    int num_shared_ways = 0;
    int a = 0; // ith node
    int b = 0; // ith + 1 node

    // Iterate through the <node_ids>
    for (int i = 0; i < size - 1; i++) 
    {
        // Initialize nodes a and b whose travel time will be found (a->b)
        a = node_ids[i];
        b = node_ids[i + 1];
        struct node * node_a = find_node(m, a);
        struct node * node_b = find_node(m, b);

        //Condition 1: Check if nodes are valid
        if (node_a == NULL) 
        { 
            printf ("error: node %d does not exist.\n", a);
            return -1.0;
        } 
        else if (node_b == NULL) 
        {
            printf("error: node %d does not exist.\n", b);
            return -1.0;
        }
        
        //Condition 2: Make sure a shared way/road exists
        num_shared_ways = node_a->num_ways;
        int * shared_way_ids = find_duplicates(node_a->way_ids, node_b->way_ids, &num_shared_ways, node_b->num_ways);
        
        if (num_shared_ways == 0) 
        {
            printf("error: there are no roads between node %d and node %d.\n", a, b);
            return -1.0;
        }

        struct way * shared_ways [num_shared_ways];
        
        // Find the ways that contain both nodes a and b
        for (int i = 0; i < num_shared_ways; i++) 
        {
            struct way * new_way = find_way(m, shared_way_ids[i]);
            shared_ways[i] = new_way;
        }

        // Free allocated memory
        free(shared_way_ids); 

        bool oneway = false;
        float maxspeed = 0.0;
        bool adjacent = false;
        struct way * curr_way;
        int order[2] = {-1,-1};
        int * nodes;

        //Condition 3: Make sure nodes are adjacent in a given array
        for (int i = 0; i < num_shared_ways; i++) 
        {
            curr_way = shared_ways[i];
            
            // Check if nodes are adjacent by iterating through the nodes of curr_way
            for (int j = 0; j < curr_way->num_nodes - 1; j++) 
            {
                nodes = curr_way->node_ids;
                
                if ((nodes[j] == a && nodes[j + 1] == b) || (nodes[j] == b && nodes[j + 1] == a)) 
                {
                    maxspeed = curr_way->maxspeed;
                    oneway = curr_way->oneway;
                    
                    // Save the order of the nodes (used to check for oneway)
                    if (nodes[j] == a) 
                    {
                        order[0] = a;
                        order[1] = b;
                    }
                    else 
                    {
                        order[0] = b;
                        order[1] = a;
                    }

                    adjacent = true;
                    break;
                }
            }
        }

        if (!adjacent) 
        {
            printf("error: cannot go directly from node %d to node %d.\n", a, b);
            return -1.0;
        }

        //Condition 4: Check if the way is oneway
        if (oneway) 
        {
            // Use the order we previously got from checking adjacent nodes to
            // check if the order is correct
            if (order[0] == b) 
            {
                printf("error: cannot go in reverse from node %d to node %d.\n", a, b);
                return -1.0;
            }
        }
                    
        //Condition 5: Check for same node
        for (int j = i + 1; j < size; j++) 
        {
            if (a == node_ids[j]) 
            {
                printf("error: node %d appeared more than once.\n", a);
                return -1.0;
            }
        }

        total_travel_time +=  distance_between_nodes(node_a, node_b) / (maxspeed / 60.0);         
    }
    
    return total_travel_time;
}


/**
 * Creates a new vertex
 * 
 * @param node The node which the vertex will be associated with 
 * @param dist The distance from the start node
 * @param prev The previous node
 * @return An intialized vertex corresponding to <n>
 */ 
struct vertex * 
create_vertex(struct node * node, double dist, struct vertex *prev) 
{
    struct vertex * new_vertex = malloc(sizeof(struct vertex));
    
    if (new_vertex == NULL) 
    {
        free(new_vertex);
        return NULL;
    }

    if (new_vertex != NULL) 
    {
        new_vertex->node = node;
        new_vertex->dist = dist;
        new_vertex->prev = prev;
    }

    return new_vertex;
}


/**
 * Restores the min-heap property to a given heap starting at a given index
 * 
 * @param heap The heap to have min_heapify performed on
 * @param size The size of the heap
 * @param index The index to start min_heapify
 */ 
void 
min_heapify(struct vertex ** heap, int size, int index) 
{
    int min = index;
    int left_child = 2 * index + 1;
    int right_child = 2 * index + 2;

    // Check if the children of heap[min] are smaller than it
    if (left_child < size && heap[left_child]->dist < heap[min]->dist) 
    {
        min = left_child;
    }

    if (right_child < size && heap[right_child]->dist < heap[min]->dist) 
    {
        min = right_child;
    }

    // If the min element is not the current node, swap the min element and call min_heapify again
    if (min != index) 
    {
        struct vertex * temp = heap[index];
        heap[index] = heap[min];
        heap[min] = temp;
        min_heapify(heap, size, min);
    }
}


/**
 * Restores the min-heap property to a given heap starting at a given index
 * 
 * @param heap The heap to extract from
 * @param size The size of the heap
 */ 
void 
remove_min_vertex(struct vertex ** heap, int * size) 
{
    (* size)--;

    // Replace root with last element
    heap[0] = heap[* size];
    
    // Reheapify
    min_heapify(heap, * size, 0);
}


// Creates and prints the shortest path from <start_id> to <end_id>
void 
ssmap_path_create(const struct ssmap * m, int start_id, int end_id) 
{
    // Check if start and end are valid
    struct node * start_node = find_node(m, start_id);
    struct node * end_node = find_node(m, end_id);
    
    if (start_node == NULL) 
    {
        printf("error: node %d does not exist.\n", start_id);
        return;
    } 
    else if (end_node == NULL) 
    {
        printf("error: node %d does not exist.\n", end_id);
        return;
    }

    // Initialize verticies to find the path
    struct vertex * verticies[m->num_nodes];
    
    for (int i = 0; i < m->num_nodes; i++) 
    {
        if (m->nodes[i]->id == start_id) 
        {
            verticies[i] = create_vertex(m->nodes[i], 0, NULL);
        }
        else
        {
            verticies[i] = create_vertex(m->nodes[i], MAX, NULL); 
        }
    }

    // Initialize the priority queue
    struct vertex * priority_q[m->num_nodes];
    priority_q[0] = create_vertex(start_node, 0, NULL);
    int pq_size = 1;

    // Dijkstra's algorithm
    while (pq_size > 0) 
    {
        // Extract min vertex
        struct vertex * curr_vertex = priority_q[0];
        remove_min_vertex(priority_q, &pq_size);
        min_heapify(priority_q, pq_size, 0);
        curr_vertex->visited = true;

        // Look through all associated ways
        for (int i = 0; i < curr_vertex->node->num_ways; i++) 
        {
            struct way * curr_way = find_way(m, curr_vertex->node->way_ids[i]);
            
            // Forward traversal
            for (int j = 0; j < curr_way->num_nodes - 1; j++) 
            {
                if (curr_way->node_ids[j] == curr_vertex->node->id) 
                {
                    // Find next node
                    struct node * curr_node = find_node(m, curr_way->node_ids[j + 1]);
                    
                    // Check if visited
                    if (curr_node != NULL && !verticies[curr_node->id]->visited) 
                    { 
                        // Calculate total distance
                        double total_dist = curr_vertex->dist + distance_between_nodes(curr_vertex->node, curr_node);
                        
                        // Check if distance is smaller if so replace
                        if (total_dist < verticies[curr_node->id]->dist) 
                        {
                            verticies[curr_node->id]->dist = total_dist;
                            verticies[curr_node->id]->prev = curr_vertex;
                            bool exists_in_q = false;
                            
                            // Check if node already exists in priority queue
                            for (int k = 0; k < pq_size; k++) 
                            {
                                if (priority_q[k]->node->id == curr_node->id) 
                                {
                                    exists_in_q = true;
                                    break;
                                }
                            }
                            
                            // If it doesn't exist in the q add it
                            if (!exists_in_q) 
                            {
                                priority_q[pq_size++] = create_vertex(curr_node, total_dist, curr_vertex);
                                min_heapify(priority_q, pq_size, pq_size - 1);
                            }
                        }
                    }
                }
            }

            // Backwards traversal
            for (int j = curr_way->num_nodes - 1; j > 0; j--) {
                if (curr_way->node_ids[j] == curr_vertex->node->id) 
                {
                    // Find previous node
                    struct node * curr_node = find_node(m, curr_way->node_ids[j - 1]);
                    
                    // Make sure way is not one-way and not visited
                    if (curr_node != NULL && !verticies[curr_node->id]->visited && (!curr_way->oneway || j == curr_way->num_nodes - 1)) 
                    { 
                        // Calculate total distance
                        double total_dist = curr_vertex->dist + distance_between_nodes(curr_vertex->node, curr_node);
                        
                        // Check if distance is smaller if so replace
                        if (total_dist < verticies[curr_node->id]->dist) 
                        {
                            verticies[curr_node->id]->dist = total_dist;
                            verticies[curr_node->id]->prev = curr_vertex;
                            bool exists_in_q = false;
                            
                            // Check if node already exists in priority queue
                            for (int k = 0; k < pq_size; k++)
                            {
                                if (priority_q[k]->node->id == curr_node->id) 
                                {
                                    exists_in_q = true;
                                    break;
                                }
                            }

                            // If it doesn't exist in the q add it
                            if (!exists_in_q) 
                            {
                                priority_q[pq_size++] = create_vertex(curr_node, total_dist, curr_vertex);
                                min_heapify(priority_q, pq_size, pq_size - 1);
                            }
                        }
                    }
                }
            }
        }
    }

    int path[m->num_nodes];
    int path_count = 0;
    
    // Print results
    if (verticies[end_id]->dist != MAX) 
    {
        struct vertex * temp = verticies[end_id];
        
        while (temp != NULL) 
        {
            path[path_count] = temp->node->id;
            temp = temp->prev;
            path_count++;
        }

        for (int i = path_count - 1; i > 0; i--) 
        {
            printf("%d ", path[i]);
        }

        printf("%d\n", path[0]);
    } 
    else 
    { 
        printf("error: could not find a path from node %d to node %d.", start_id, end_id);
    }

    // Free everything
    for (int i = 0; i < pq_size; i++) 
    {
        free(priority_q[i]);
    }
}