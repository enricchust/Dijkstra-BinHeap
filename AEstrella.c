// 1569251 Sergi Cant�n Sim�

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>
# include <limits.h>

# define R 6371000
# define TORAD 3.1415926536/180
# define MAX_NUM_ARESTES 10
# define tipusIdNode long int
# define INFINIT 1.79769e+308

typedef struct {
    char idCarrer[12];
    int numSeguentNode;
} Aresta;

typedef struct ElementLlista {
    tipusIdNode id, idPare;
    double latitud, longitud, cost;
    int numArestes;
    Aresta arestes[MAX_NUM_ARESTES];
    struct ElementLlista *seguent;
} Node;


double distancia(Node*, Node*);
int buscaNode(Node[], int, tipusIdNode);
void insereixNodeLlista(Node**, Node*);
void esborraNodeLlista(Node**, tipusIdNode);
int nodeEsALlista(Node**, long long int);
tipusIdNode buscaNodeMenorCost(Node**, Node*);
tipusIdNode* AEstrella(tipusIdNode, tipusIdNode, Node[], int);


int main(int argc, char *argv[]) {
    FILE *fitxerNodes;
    FILE *fitxerCarrers;
    int numNodes = 0, numCarrers = 0;
    char c;
    Node *nodes; // Llista on es guarda la informaci� de cada node
    tipusIdNode idNodeInicial, idNodeFinal;

    if (argc != 3) {
        printf("\nERROR. El nombre d'arguments es diferent de 2.\n");
        return 3;
    } else {
        idNodeInicial = atoll(argv[1]);
        idNodeFinal = atoll(argv[2]);
    }

    // Lectura del fitxer Nodes.csv ---------------
    fitxerNodes = fopen("Nodes.csv", "r");
    if (fitxerNodes == NULL) {
        printf("\nERROR. No s'ha pogut accedir al fitxer Nodes.csv\n");
        return 1;
    }
    while ((c = fgetc(fitxerNodes)) != EOF) {
        if(c == ';')  numNodes++;
    }
    numNodes /= 2;
    if((nodes = (Node*) malloc(numNodes*sizeof(Node))) == NULL) {
        printf("ERROR. No s'ha pogut assignar suficient memoria.");
    }
    rewind(fitxerNodes);
    for (int i = 0; i < numNodes; i++) {
        fscanf(fitxerNodes, "%ld;%lf;%lf", &(nodes[i].id), &(nodes[i].latitud), &(nodes[i].longitud));
        nodes[i].numArestes = 0;
        nodes[i].cost = INFINIT;
    }
    fclose(fitxerNodes);
    //---------------------------------------------

    if(buscaNode(nodes, numNodes, idNodeInicial) == -1 || buscaNode(nodes, numNodes, idNodeFinal) == -1) {
        printf("ERROR. No s'han trobat els nodes amb les id introduides.\n");
        return 4;
    }

    // Lectura del fitxer Carrers.csv -------------
    fitxerCarrers = fopen("Carrers.csv", "r");
    if (fitxerCarrers == NULL) {
        printf("\nNo s'ha accedit al fitxer Carrers.csv\n");
        return 2;
    }
    while ((c = fgetc(fitxerCarrers))!=EOF) {
        if(c == '=')  numCarrers++;
    }
    int posicioNode1, posicioNode2;
    tipusIdNode idNode;
    char idCarrer[12];

    rewind(fitxerCarrers);

    for (int i = 0; i < numCarrers; i++) {
        while ((c = fgetc(fitxerCarrers)) != '=');
        fscanf(fitxerCarrers, "%10s", idCarrer);
        fgetc(fitxerCarrers);
        fscanf(fitxerCarrers, "%ld", &idNode);
        posicioNode1 = buscaNode(nodes, numNodes, idNode);

        do {
            fgetc(fitxerCarrers); //treu els ;
            idNode = -1;
            fscanf(fitxerCarrers, "%ld", &idNode);

            if(idNode != -1) {
                posicioNode2 = buscaNode(nodes, numNodes, idNode);
                strcpy(nodes[posicioNode1].arestes[nodes[posicioNode1].numArestes].idCarrer, idCarrer);
                nodes[posicioNode1].arestes[nodes[posicioNode1].numArestes].numSeguentNode = posicioNode2;
                nodes[posicioNode1].numArestes++;
                strcpy(nodes[posicioNode2].arestes[nodes[posicioNode2].numArestes].idCarrer, idCarrer);
                nodes[posicioNode2].arestes[nodes[posicioNode2].numArestes].numSeguentNode = posicioNode1;
                nodes[posicioNode2].numArestes++;
                posicioNode1 = posicioNode2;
            }
        } while (idNode != -1);
    }
    fclose(fitxerCarrers);
    //---------------------------------------------

    tipusIdNode *cami = AEstrella(idNodeInicial, idNodeFinal, nodes, numNodes);
    int numNodesCami = 0;
    while (cami[numNodesCami] != -1) {
        numNodesCami++;
    }
    double distanciaOptima = nodes[buscaNode(nodes, numNodes, cami[numNodesCami-1])].cost;
    printf("\nLa distancia de %s a %s es de %lf metres.\n", argv[1], argv[2], distanciaOptima);
    printf("Cami optim:\n");
    printf("-----------------------------------------------------------------\n");
    Node *nodeCami;
    for (int i=0; i < numNodesCami; i++) {
        nodeCami = &nodes[buscaNode(nodes, numNodes, cami[i])];
        printf("Id=");
        if(nodeCami->id < 10e8) {
            printf("0");
        }
        printf("%ld | %lf | %lf | Dist=%lf\n", nodeCami->id, nodeCami->latitud, nodeCami->longitud, nodeCami->cost);
    }
    printf("-----------------------------------------------------------------\n");
}


double distancia(Node *node1, Node *node2){
    double x1 = R * cos(node1->longitud*TORAD)*cos(node1->latitud*TORAD);
    double y1 = R * sin(node1->longitud*TORAD)*cos(node1->latitud*TORAD);
    double z1 = R * sin(node1->latitud*TORAD);
    double x2 = R * cos(node2->longitud*TORAD)*cos(node2->latitud*TORAD);
    double y2 = R * sin(node2->longitud*TORAD)*cos(node2->latitud*TORAD);
    double z2 = R * sin(node2->latitud*TORAD);
    double distance = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));
    return distance;
}


int buscaNode(Node nodes[], int numNodes, tipusIdNode id){
        unsigned primer = 0, ultim = numNodes;
        unsigned mitja = (primer + ultim) / 2;
        unsigned i;

        while(primer <= ultim){
            if(nodes[mitja].id < id){
                primer = mitja + 1;
            } else if(nodes[mitja].id == id){
                i = mitja;
                break;
            }else {
                ultim = mitja - 1;
            }
            mitja = (primer + ultim) / 2;
        }
        return i;
    }
void insereixNodeLlista(Node **iniciLlista, Node *actual) {
    if(*iniciLlista == NULL) {
        *iniciLlista = actual;
        actual->seguent = NULL;
    } else {
            actual->seguent = *iniciLlista;
            *iniciLlista = actual;
    }
}

void esborraNodeLlista(Node **iniciLlista, tipusIdNode idNodeAEsborrar) {
    Node *anterior = *iniciLlista;
    if (anterior->id == idNodeAEsborrar) {
        *iniciLlista = anterior->seguent;
    } else {
        while(anterior->seguent->id != idNodeAEsborrar) {
            if (anterior->seguent == NULL) {
                break; // No s'ha trobat el node
            }
            anterior=anterior->seguent;
        }
        anterior->seguent = anterior->seguent->seguent;
    }
}


int nodeEsALlista(Node **iniciLlista, long long int id) {
    Node *actual = *iniciLlista;
    if(*iniciLlista!=NULL) {
        while(actual != NULL){
            if(actual->id == id){
                return 1;
            }
            actual = actual->seguent;
        }
    }
    return 0; // El node no �s a la llista
}

tipusIdNode buscaNodeMenorCost(Node **iniciLlista, Node *nodeFinal) {
    Node *actual = *iniciLlista;
    double distanciaHeuristica = distancia(actual, nodeFinal);
    tipusIdNode idNodeMenorCost = -1;
    double menorCost = actual->cost + distanciaHeuristica;
    idNodeMenorCost = actual->id;
    while(actual != NULL) {
        if(actual->cost + distanciaHeuristica < menorCost) {
            menorCost = actual->cost + distanciaHeuristica;
            idNodeMenorCost = actual->id;
        }
        actual = actual->seguent;
    }
    return idNodeMenorCost;
}

tipusIdNode* AEstrella(tipusIdNode idNodeInicial, tipusIdNode idNodeFinal, Node nodes[], int numNodes) {
    Node *llistaOberta = NULL;
    Node *llistaTancada = NULL;
    Node *nodeInicial = &nodes[buscaNode(nodes, numNodes, idNodeInicial)];
    Node *nodeFinal = &nodes[buscaNode(nodes, numNodes, idNodeFinal)];
    Node *nodeActual;
    Node *nodeSuccessor;
    double costActualSuccessor;
    int numNodeSuccessor;
    tipusIdNode idNodeMenorCost;
    nodeInicial->cost = 0;
    insereixNodeLlista(&llistaOberta, nodeInicial);

    while (llistaOberta != NULL) {
        idNodeMenorCost = buscaNodeMenorCost(&llistaOberta, nodeFinal);
        nodeActual = &nodes[buscaNode(nodes, numNodes, idNodeMenorCost)];
        esborraNodeLlista(&llistaOberta, idNodeMenorCost);

        if (nodeActual->id == idNodeFinal) break; // Soluci� trobada

        int numExpansions = nodeActual->numArestes;
        for (int i = 0; i < numExpansions; i++) {
            numNodeSuccessor = nodeActual->arestes[i].numSeguentNode;
            nodeSuccessor = &nodes[numNodeSuccessor];
            costActualSuccessor = nodeActual->cost + distancia(nodeActual, nodeSuccessor);
            if (nodeEsALlista(&llistaOberta, nodeSuccessor->id)) {
                if (nodeSuccessor->cost <= costActualSuccessor) continue;
            } else if (nodeEsALlista(&llistaTancada, nodeSuccessor->id)) {
                if (nodeSuccessor->cost <= costActualSuccessor) continue;
                esborraNodeLlista(&llistaTancada, nodes[numNodeSuccessor].id);
                insereixNodeLlista(&llistaOberta, nodeSuccessor);
            } else {
                insereixNodeLlista(&llistaOberta, nodeSuccessor);
            }
            nodeSuccessor->cost = costActualSuccessor;
            nodeSuccessor->idPare = nodeActual->id;
        }
        insereixNodeLlista(&llistaTancada, nodeActual);
    }

    if (nodeActual->id != nodeFinal->id) {
        printf("ERROR. No s'ha trobat cap cami.\n");
    } else {
        // Ja a'ha arribat al node final. Ara cal refer el cam� al rev�s per guardar i retornar la informaci� dels nodes que el formen

        int numNodesCami = 1;
        Node *auxiliar = nodeFinal;
        auxiliar = &nodes[buscaNode(nodes, numNodes, auxiliar->idPare)];
        while(auxiliar->id != nodeInicial->id){
            auxiliar = &nodes[buscaNode(nodes, numNodes, auxiliar->idPare)];
            numNodesCami++;
        }

        tipusIdNode *idNodesCami;
        if ((idNodesCami = (tipusIdNode *) malloc(sizeof(tipusIdNode) * (numNodesCami + 1))) == NULL) {
            printf("ERROR. No s'ha pogut assignar suficient memoria.\n");
            exit(1);
        }
        auxiliar = nodeFinal;

        for (int i = 0; i < numNodesCami; i++) {
            idNodesCami[numNodesCami-i-1] = auxiliar->id;
            if (i != numNodesCami - 1) {
                auxiliar = &nodes[buscaNode(nodes, numNodes, auxiliar->idPare)];
            }
        }
        idNodesCami[numNodesCami] = -1;
        return idNodesCami;
    }
}

