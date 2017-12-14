//
//  fire.h
//  cultimate
//
//  Created by Jiahao Zhang on 12/12/17.
//  Copyright © 2017 Jiahao Zhang. All rights reserved.
//

#ifndef fire_h
#define fire_h
#include "particle.h"
void fire(int,double,particle*);
void addptolist(parnode* head,parnode* tail,int in);
void printplist(parnode* head);
void updatepincell(int N,double len,double celllen,particle* allpart);
void updatecellofp(int N,int cellsize,parnode* cellall[],particle *allpart);
double distance(int index1,int index2,double len,particle* allpart);
void updateoneneigh(int ind,int cellsize,double len,double rverlet,parnode* cellall[],particle* allpart);
void updateallneigh(int N,int cellsize,double len,double rverlet,parnode* cellall[],particle*);
void updateforce(int N,double len,particle* allpart);
double* leapfrogone(int N,double Dt,double len,particle* allpart);
void leapfrogtwo(int N,double Dt,particle* allpart);
void setv(int N,double alpha,particle* allpart);
double power(int N,particle* allpart);
void freeze(int N,particle* allpart);
double energy(int N,double len,particle* allpart);
double kinetic(int N,particle* allpart);
int* bias(int index,int cellsize,particle* allpart);
#endif /* fire_h */
