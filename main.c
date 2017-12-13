//
//  main.c
//  cultimate
//
//  Created by Jiahao Zhang on 12/12/17.
//  Copyright © 2017 Jiahao Zhang. All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fire.h"
#include "myrandom.h"
#include "particle.h"
#include <time.h>
#define pi 3.1415926535897932385
int main(int argc, const char * argv[]) {
    time_t t1,t2;
    time(&t1);
    double fraction=0.8;
    int N=4000;
    double volume=0;
    double temp=0;
    double len;
    particle allpart[N];
    //***********start initial system for particles***************//
    for (int i=0; i<N; i++) {
        allpart[i].radius=1;
        allpart[i].neighbor=(parnode*)malloc(sizeof(parnode));
        allpart[i].tail=allpart[i].neighbor;
        allpart[i].neighbor->index=-1;
				allpart[i].neighbor->next=NULL;
    }
    for (size_t i=0; i<N; i++) {
        temp=allpart[i].radius;
        volume=volume+pi/3.0*4*temp*temp*temp;
    }
    len=cbrt(volume/fraction);
    for (size_t i=0; i<N; i++) {
        for (size_t j=0; j<3; j++) {
            allpart[i].posit[j]=len*(genrand()-0.5);
            allpart[i].speed[j]=0.0;
            allpart[i].force[j]=0.0;
        }
    }
    //**************end initialization for particles**************//
    fire(N, len, allpart);
    time(&t2);
    printf("the length of the system is: %lf time used: %lf\n",len,difftime(t2, t1));
    return 0;
}
