//Enric Chust Gimeno 1638322

#include <stdio.h>
#include <stdlib.h>
#include <values.h> // Per MAXFLOAT
#include <math.h>

void ExitError(const char *lloc, const char *miss, int err) { fprintf (stderr, "\nERROR %s: %s.\nAbortem el procés ...\n\n", lloc, miss); exit(err); }

typedef struct{ unsigned num_node_adjacent; float temps; } aresta;
typedef struct{ unsigned n_arestes; aresta arestes[5]; } node;
typedef struct{ float cost_g; unsigned node_pare; } node_routing_status;

#define NNODES 21

int Dijkstra(unsigned, node *, node_routing_status *, unsigned);

int main() { register unsigned v;
/* Inicialitzacio */
    node_routing_status status_nodes[NNODES];
    node nodes[NNODES] = {
            {3, { {1, 0.528},   {2, 0.495},   {3, 0.471} }},                 // 0
            {2, { {0, 0.528},   {3, 0.508} }},                               // 1
            {4, { {0, 0.495},   {3, 3.437},   {5, 12.033},  {15, 34.852} }}, // 2
            {4, { {0, 0.471},   {1, 0.508},   {2, 3.437},   {4, 23.155} }},  // 3
            {4, { {3, 23.155},  {5, 6.891},   {6, 4.285},   {7, 0.520} }},   // 4
            {2, { {2, 12.033},  {4, 6.8910} }},                              // 5
            {3, { {4, 4.285},   {7, 0.630},   {8, 17.406} }},                // 6
            {2, { {4, 0.520},   {6, 0.630} }},                               // 7
            {5, { {6, 17.406},  {9, 6.657},   {10, 15.216}, {11, 10.625}, {12, 17.320} }}, // 8
            {2, { {8, 6.657},   {12, 16.450} }},                             // 9
            {2, { {8, 15.216},  {14, 12.373} }},                             // 10
            {2, { {8, 10.625},  {12, 3.618} }},                              // 11
            {3, { {8, 17.320},  {9, 16.450},  {11, 3.618} }},                // 12
            {2, { {14, 4.450},  {19, 6.450} }},                              // 13
            {4, { {10, 12.373}, {13, 4.450},  {15, 16.178}, {19, 5.203} }},  // 14
            {5, { {2, 34.852},  {14, 16.178}, {16, 4.818},  {18, 3.877},  {19, 19.131} }}, // 15
            {3, { {15, 4.818},  {17, 3.199},  {18, 2.976} }},                // 16
            {2, { {16, 3.199},  {18, 20.832} }},                             // 17
            {4, { {15, 3.877},  {16, 2.976},  {17, 20.832}, {20, 2.510} }},  // 18
            {4, { {13, 6.450},  {14, 5.203},  {15, 19.131}, {20, 13.313} }}, // 19
            {2, { {18, 2.510},  {19, 13.313} }} };                           // 20
    unsigned node_origen = 0;

    if(Dijkstra(node_origen, nodes, status_nodes, NNODES) == 66) ExitError("a Dijkstra", "no és possible allocatar memòria per la cua", 66);

/* Detecció dels nodes extremals: són els que no són pare de ningú */
    char EsPare[NNODES] = { 0U };
    for (v=1; v < NNODES; v++) {
        if(!(status_nodes[v].cost_g < MAXFLOAT)) ExitError("al graf", "hi ha nodes inaccessibles", 3);
        EsPare[status_nodes[v].node_pare] = 1;
    }

/* Per cada node extremal trobem el camí des-de l'origen */
    unsigned Cami[NNODES];
    for(v=1; v < NNODES; v++) {
        if(EsPare[v]) continue;
        register int l = 0;
        Cami[0] = v;
        do { l++; Cami[l] = status_nodes[Cami[l-1]].node_pare;} while (Cami[l]) ;
        fprintf(stdout, "[0]");
        for(l-- ; l >= 0; l--) fprintf(stdout, " -- %g --> [%d]", status_nodes[Cami[l]].cost_g, Cami[l]) ;
        fprintf(stdout, "\n");
    }
    return 0;
}

#define MAXNumlevels 60
typedef struct {
    short nlevels;
    unsigned long last_level_size;
    unsigned *Cua_Prio_BiHe_level[MAXNumlevels];
} Cua_Prio_BiHe;

typedef struct {
    int level;
    unsigned long index;
} CPBH_Element;

int CPBH_EsBuit(Cua_Prio_BiHe *);
unsigned CPBH_desencua(Cua_Prio_BiHe *, node_routing_status *);
int CPBH_encua(unsigned, node_routing_status *, Cua_Prio_BiHe *);
void CPBH_reencua(unsigned, node_routing_status *, Cua_Prio_BiHe *);

void CPBH_heapify_up(CPBH_Element node, node_routing_status *nods, Cua_Prio_BiHe *BH);
void CPBH_heapify_down(CPBH_Element e, node_routing_status *nods, Cua_Prio_BiHe *BH);
double CPBH_getcost(CPBH_Element node, node_routing_status *nods, Cua_Prio_BiHe *BH);
CPBH_Element CPBH_pare(CPBH_Element node, Cua_Prio_BiHe *BH);
CPBH_Element CPBH_right_son(CPBH_Element node, Cua_Prio_BiHe *BH);
CPBH_Element CPBH_left_son(CPBH_Element node, Cua_Prio_BiHe *BH);


/* La fase d'inicialització de Dijkstra inclou establir
 *          nodstatus[u].cost_g := MAXFLOAT;     i     Node_es_a_S[v] = 0;
 * per a tot node v.
 *
 * Per altra banda quan un node entra a la cua (s'encua) ho fa amb nodstatus[v].cost_g < MAXFLOAT
 * i pot fer una de les dues accions següents:
 *     * Moure's dins de la cua (reencua) amb un valor nodstatus[v].cost_g < MAXFLOAT més baix del que tenia fins la iteració actual; o
 *     * Sortir de la cua amb la regla ExtractMin (Desencua). En aquest cas sabem que el seu cost ja és optimal i l'hem d'afegir al conjunt S
 *       En aquest cas tenim nodstatus[u].cost_g < MAXFLOAT i marquem que v ∈ S establint Node_es_a_S[v] = 1;
 *
 * Llavors, un node arbitrari v pot estar en un dels tres estats següents
 * v ∈ S (S és el conjunt de nodes dels que ja se sap la distància mínima) <==> Node_es_a_S[v] == 1 (notis que, en aquest cas, tenim nodstatus[u].cost_g < MAXFLOAT)
 * v ∉ S  i, a més, v no és a la Cua_Prio_BiHe <==> nodstatus[u].cost_g == MAXFLOAT (notis que, en aquest cas, tenim Node_es_a_S[v] == 0)
 * v ∉ S però v és a la Cua_Prio_BiHe <==> nodstatus[u].cost_g < MAXFLOAT && Node_es_a_S[v] == 0 */

int Dijkstra(unsigned source, node *nods, node_routing_status *nodstatus, unsigned nnods ){ register unsigned v;
    char Node_es_a_S[nnods];
    for (v=0; v < nnods; v++) { nodstatus[v].cost_g = MAXFLOAT; Node_es_a_S[v] = 0; } // Totes les distàncies a infinit

    Cua_Prio_BiHe LaCua = { 0, 0 };
    nodstatus[source].cost_g = 0.0; nodstatus[source].node_pare = UINT_MAX; // El node source és l'arrel i, per tant, no té pare
    CPBH_encua(source, nodstatus, &LaCua);

// Llaç principal
    while(!CPBH_EsBuit(&LaCua)){ register unsigned e;
        unsigned u = CPBH_desencua(&LaCua, nodstatus); Node_es_a_S[u] = 1;
        for(e=0; e < nods[u].n_arestes; e++){ /* Llaç d'expansió del node u que acabem de posar a S */
            v = nods[u].arestes[e].num_node_adjacent;
            if(Node_es_a_S[v]) continue;
            float cost = nodstatus[u].cost_g + nods[u].arestes[e].temps;
            if(nodstatus[v].cost_g > cost) {
                int IsNodeInLaCua = nodstatus[v].cost_g < MAXFLOAT; /* Relaxation step */
                nodstatus[v].cost_g = cost; nodstatus[v].node_pare = u;
                if(IsNodeInLaCua) CPBH_reencua(v, nodstatus, &LaCua); else { int res_encua = CPBH_encua(v, nodstatus, &LaCua); if(res_encua) return 66; }
            } /* Fi del Relaxation step */
        } /* Fi del llaç d'expansió del node u que acabem de posar a S */
    } /* Fi del Llaç principal */
    return 0;
}

/***********************************
 * Funcions de cua amb Binary Heap *
 ********************************************************************************
 * Recordem la Propietat de forma del BinHeap:
 * el BinHeap és un arbre binari complet; és a dir, tots els nivells de l’arbre,
 * excepte possiblement l’últim (el més profund) estan completament omplerts
 * (amb dos fills per node) i, si l’últim nivell de l’arbre no està complet,
 * els nodes d’aquest nivell s’omplen d’esquerra a dreta.
 *
 * La Propietat de forma del BinHeap traduida a la representació
 *      typedef struct {
 *          char nlevels;
 *          unsigned long last_level_size;
 *          unsigned *Cua_Prio_BiHe_level[MAXNumlevels];
 *      } Cua_Prio_BiHe;
 * s'implmenta d'acord amb la següent
 * CONVENCIÓ: Tots els nivells l (de 0 a nlevels-2) han d'estar plens (amb 2^l elements).
 *            L'últim nivell (l = nlevels-1) ha de ser no buit, però pot no estar complet ==>
 *                    el nombre d'elements del darrer nivell (nlevels-1), denotat per last_level_size,
 *                    ha de verificar 0 < last_level_size <= 2^(nlevels-1)
 *            nlevels = 0 és un heap sense nivells; és a dir, un heap buit
 *            OBSERVEU que, quan el heap és no buit, el nivell 0 sempre està ple, ja que té com a màxim i com a mínim un element.
* MAXNumlevels = 60: Del que s'ha dit abans, tenim que la capacitat màxima total del heap és de 2^(MAXNumlevels)-1 = 2^(60)-1, realment més que suficient per a una cua.
 *
 * A més cal que es compleixi la Propietat Heap:
 * la clau emmagatzemada a cada node és menor o igual que les claus dels fills del node, segons algun ordre total.
 * Iterativament la Propietat heap diu que el cost d'un node donat N és més petit o igual que els valors de tots els nodes del subarbre que té el node N com arrel.
 ******************************************************************************************************************************************************************/


int CPBH_EsBuit(Cua_Prio_BiHe *BH){return(BH->nlevels == 0);}

int CPBH_encua(unsigned numNode, node_routing_status *nods, Cua_Prio_BiHe *BH){
    if(CPBH_EsBuit(BH)){
        //si es buit creem el primer nivell i li assignem la memoria
        BH->nlevels = 1;
        BH->last_level_size = 1;
        BH->Cua_Prio_BiHe_level[0] = (unsigned *)malloc((pow(2,0))*sizeof(unsigned));
        if (BH->Cua_Prio_BiHe_level[0] == NULL){
            printf("No és possible assignar la memòria necessària.\n");
            return 1;
        }
        BH->Cua_Prio_BiHe_level[0][0] = numNode;
    }else {
        if (BH->last_level_size == pow(2,((BH->nlevels)-1))) { //mirem si l'ultim nivell ja esta ple, en cas de ser-ho, creem un nou nivell,
            BH->nlevels += 1;
            BH->last_level_size = 1;
            BH->Cua_Prio_BiHe_level[BH->nlevels-1] = (unsigned *)malloc((pow(2,(BH->nlevels-1)))*sizeof(unsigned));
            if (BH->Cua_Prio_BiHe_level[BH->nlevels-1] == NULL){
                printf("No és possible assignar la memòria necessària.\n");
                return 1;
            }
            BH->Cua_Prio_BiHe_level[BH->nlevels-1][BH->last_level_size-1] = numNode;
        } else { //si no esta ple indiquem que l'ultim nivell tindra un nou node
            BH->last_level_size += 1;
            BH->Cua_Prio_BiHe_level[BH->nlevels-1][BH->last_level_size-1] = numNode;
        }
        //afegim al node el seu nivell i index
        CPBH_Element node = {BH->nlevels - 1, BH->last_level_size - 1};

        CPBH_heapify_up(node, nods, BH);
    }
    return 0;
}

void CPBH_reencua(unsigned numNode, node_routing_status *nods, Cua_Prio_BiHe *BH){
    CPBH_Element node;
    //recorrem l'arbre en busca del node que hem de reencuar, un vegada l'haguem trobat fem us de la funcio heapify_up com amb encua
    for(int level=0; level<BH->nlevels; level++){
        for(unsigned long index=0; index<(pow(2,level)); index++){
            if(BH->Cua_Prio_BiHe_level[level][index] == numNode){
                node.level = level;
                node.index = index;
                break;
            }
        }
    }
    CPBH_heapify_up(node, nods, BH);
}

unsigned CPBH_desencua(Cua_Prio_BiHe *BH, node_routing_status *nods){
    if(CPBH_EsBuit(BH)){
        return 1;
    }
    unsigned v = BH->Cua_Prio_BiHe_level[0][0];

    int numNode = BH->Cua_Prio_BiHe_level[BH->nlevels - 1][BH->last_level_size - 1];
    if(BH->nlevels == 1){
        free(BH->Cua_Prio_BiHe_level[0]);
        BH->nlevels = 0;
        BH->last_level_size = 0;
        return v;
    }else {
        BH->Cua_Prio_BiHe_level[0][0] = numNode;
        if (BH->last_level_size > 1) { //mirem si hi ha mes d'un element a l'ultim nivell
            BH->last_level_size -= 1;
        } else { //si nomes hi ha un element, haurem d'eliminar tambe el nivell
            BH->nlevels -= 1;
            BH->last_level_size = pow(2, (BH->nlevels - 1));
            free(BH->Cua_Prio_BiHe_level[BH->nlevels]);
        }
        CPBH_Element node = {0, 0};
        CPBH_heapify_down(node,  nods, BH);
        return v;
    }
}

void CPBH_heapify_up(CPBH_Element node, node_routing_status *nods, Cua_Prio_BiHe *BH){
    //calculem el cost del node i el node fill
    double cost = CPBH_getcost(node, nods, BH);
    CPBH_Element nodePare = CPBH_pare(node, BH);
    double costPare = CPBH_getcost(nodePare, nods, BH);

    //fem un bucle mentre no arribem al node arrel i el cost sigui menor al del pare
    while(node.level > 0 && cost < costPare){
        //intercanviem els nodes
        int numNode = BH->Cua_Prio_BiHe_level[node.level][node.index];
        BH->Cua_Prio_BiHe_level[node.level][node.index]=BH->Cua_Prio_BiHe_level[nodePare.level][nodePare.index];
        BH->Cua_Prio_BiHe_level[nodePare.level][nodePare.index]=numNode;

        //ara el node passara a tenir el nivell i index del node pare, mentre que el nou node pare sera el node pare del pare que ja es l'actual node
        node = nodePare;
        nodePare = CPBH_pare(node, BH);
        cost = CPBH_getcost(node, nods, BH);
        costPare = CPBH_getcost(nodePare, nods, BH);
    }
}

CPBH_Element CPBH_left_son(CPBH_Element node, Cua_Prio_BiHe *BH){
    CPBH_Element leftSon;
    if(node.level <= (BH->nlevels - 2)){ //ens assegurem que no sigui a l'ultim nivell i per tant que tingui fill
        leftSon.level = node.level + 1; //el nivell del fill sempre sera un mes que el del pare
        leftSon.index = node.index * 2; //l'index del fill esquerre sera parell, aixi que sera l'index del pare multiplicat per dos
    }else{ //si no te fills li assignem que el nivell es -1
        leftSon.level = -1;
    }
    return leftSon;
}

CPBH_Element CPBH_right_son(CPBH_Element node, Cua_Prio_BiHe *BH){
    CPBH_Element rightSon;
    if (node.level <= (BH->nlevels - 2)){ //ens assegurem que no sigui a l'ultim nivell i per tant que tingui fill
        rightSon.level = node.level + 1; //l'index del fill dret sera imparell, aixi que sera l'index del par multiplicat per dos mes 1
        rightSon.index = node.index * 2 + 1;
    }else{ //si no te fills li assignem que el nivell es -1
        rightSon.level = -1;
    }
    return rightSon;
}

void CPBH_heapify_down(CPBH_Element node, node_routing_status *nods, Cua_Prio_BiHe *BH){
    // calculem els diferents costos
    CPBH_Element nodeR = CPBH_right_son(node, BH);
    CPBH_Element nodeL = CPBH_left_son(node, BH);
    double cost = CPBH_getcost(node, nods, BH);
    double costL = CPBH_getcost(nodeL, nods, BH);
    double costR = CPBH_getcost(nodeR, nods, BH);

    //fem un bucle mentre no arribem al ultim nivell, ja que si arribem al ultim el node no te fills
    while(node.level < (BH->nlevels - 1)){
        if(nodeR.level == -1){
            if(costL < cost){ //comparem el cost de l'esquerra ja que no hi ha dret, si es menor s'intercanvien
                //guardem el numer de node en un altra variable ja que volem assignar el valor actual després de mmodificar-lo
                int numNode = BH->Cua_Prio_BiHe_level[node.level][node.index];
                BH->Cua_Prio_BiHe_level[node.level][node.index] = BH->Cua_Prio_BiHe_level[nodeL.level][nodeL.index];
                BH->Cua_Prio_BiHe_level[nodeL.level][nodeL.index] = numNode;
                node = nodeL; //el node que estem modificant de posicio passa a tenir el nivell i identificador del node fill
                return;
            }
        }
        else{
            if(costL < cost || costR < cost){ //si hi ha node esquerre i dret, mirem si algu dels dos te un cost menor
                if (costL <= costR){ //fem el mateix que abans segons si el cost menor es el de la dreta o el de l'esquerra
                    int numNode = BH->Cua_Prio_BiHe_level[node.level][node.index];
                    BH->Cua_Prio_BiHe_level[node.level][node.index] = BH->Cua_Prio_BiHe_level[nodeL.level][nodeL.index];
                    BH->Cua_Prio_BiHe_level[nodeL.level][nodeL.index] = numNode;
                    node = nodeL;
                }else{
                    int numNode = BH->Cua_Prio_BiHe_level[node.level][node.index];
                    BH->Cua_Prio_BiHe_level[node.level][node.index] = BH->Cua_Prio_BiHe_level[nodeR.level][nodeR.index];
                    BH->Cua_Prio_BiHe_level[nodeR.level][nodeR.index] = numNode;
                    node = nodeR;
                }
            }else{
                return;
            }
        }

        //actualitzem els nous nodes i els nous costos
        cost = CPBH_getcost(node, nods, BH);
        nodeL = CPBH_left_son(node, BH);
        nodeR = CPBH_right_son(node, BH);
        costL = CPBH_getcost(nodeL, nods, BH);
        costR = CPBH_getcost(nodeR, nods, BH);
    }
}

CPBH_Element CPBH_pare(CPBH_Element node, Cua_Prio_BiHe *BH){
    CPBH_Element nodePare;
    if(node.level > 0){
        nodePare.level = node.level - 1; //el nivell del pare sempre es un menys que el del fill
        nodePare.index = node.index / 2; //l'index del pare sempre sera la part entera del index del fill entre dos
    }
    return nodePare;
}

//primer havia intentat de fer el programa sense la funcio getcost, pero a l'hora de fer els heapify se'm complicava molt el codi
//ja que tenia que fer un us del idPare molt lios, finalment he optat per fer la funcio i calculant el cost es molt mes facil
double CPBH_getcost(CPBH_Element node, node_routing_status *nods, Cua_Prio_BiHe *BH){
    if (node.level >= 0){
        double cost = nods[BH->Cua_Prio_BiHe_level[node.level][node.index]].cost_g;
        return cost;
    }
    else{
        return MAXFLOAT;
    }
}
