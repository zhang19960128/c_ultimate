//
//  particle.h
//  cultimate
//
//  Created by Jiahao Zhang on 12/12/17.
//  Copyright © 2017 Jiahao Zhang. All rights reserved.
//

#ifndef particle_h
#define particle_h
struct parnodep{
    int index;
    struct parnodep* next;
    struct parnodep* tail;//specify the tail pointer for the list.
}parnodep;
typedef struct parnodep parnode;
struct particlep{
    double posit[3];
    double speed[3];
    double force[3];
    double cindex[3];
    double radius;
    parnode* neighbor;
    parnode* tail;//specify the tail pointer for the list
}particlep;
typedef struct particlep particle;
#endif /* particle_h */
