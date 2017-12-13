//
//  fire.c
//  cultimate
//
//  Created by Jiahao Zhang on 12/12/17.
//  Copyright © 2017 Jiahao Zhang. All rights reserved.
//
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "fire.h"
#include "particle.h"
void fire(int N,double len,particle* allpart){
    int N_min=5;
    double Dt=0.05;
    double f_inc=1.1;
    double f_dec=0.5;
    double alpha_start=0.1;
    double f_alpha=0.99;
    double celllen=0;//specify the cut off length of repusive force.
    double e_before,e_end,pow;
    double alpha=alpha_start;
    int count=0;
    double Dt_max=10*Dt;
    int cellsize=0;//specify the length of each cell.
    for (size_t i=0; i<N; i++) {
        if(celllen<allpart[i].radius) celllen=allpart[i].radius;
    }
    celllen=celllen*2;//the cutoff length should be twice the maximum radius;
    cellsize=ceil(len/celllen);//show how many cell on each edge.
    parnode *cellall[cellsize*cellsize*cellsize];
    //%%%%%%specify the head for each cell%%%%%%%%//
    for (size_t i=0; i<cellsize*cellsize*cellsize; i++) {
        cellall[i]=(parnode*)malloc(sizeof(parnode));
        cellall[i]->next=NULL;
        cellall[i]->tail=cellall[i];
        cellall[i]->index=-1;//use index=-1 to show that this is the head and useless.
    }
    //%%%%%%starting fire algorithm%%%%%%%%%%%%%%//
    int i=0;
    e_end=0;
    updatepincell(N, len, celllen, allpart);
    updatecellofp(N, cellsize, cellall, allpart);
    updateallneigh(N, cellsize, len, cellall, allpart);
    updateforce(N, len, allpart);
    do{
        i++;
        e_before=e_end;
        updatepincell(N, len, celllen, allpart);
        updatecellofp(N, cellsize, cellall, allpart);
        updateallneigh(N, cellsize, len, cellall, allpart);
        leapfrogone(N, Dt, len, allpart);
        updateforce(N, len, allpart);
        leapfrogtwo(N, Dt, allpart);
        /*
        pow=power(N,allpart);
        setv(N,alpha,allpart);
        if(pow>0&&count>N_min){
                count=0;
                Dt=Dt*f_inc<Dt_max ? Dt*f_inc:Dt_max;
                alpha=alpha*f_alpha;
        }
        else if (pow>0){
            count=0;
        }
        else{
            count++;
            Dt=Dt*f_dec;
            alpha=alpha_start;
            freeze(N,allpart);
        }
         */
        e_end=energy(N,len,allpart);
        printf("step:%d energy is: %lf\n",i,e_end+kinetic(N, allpart));
    }while(i<5000);
};
double energy(int N,double len,particle* allpart){
    parnode* temp;
    double dij;
    double rij;
    double ener=0;
    for (int i=0;i<N; i++) {
        temp=allpart[i].neighbor->next;
        while (temp!=NULL) {
            rij=distance(temp->index, i, len, allpart);
            dij=allpart[temp->index].radius+allpart[i].radius;
            ener=ener+(1-rij/dij)*(1-rij/dij);
            temp=temp->next;
        }
    }
    return ener/2;
}
double kinetic(int N,particle* allpart){
    double kall=0;
    for (int i=0; i<N; i++) {
        for (int j=0; j<3; j++) {
            kall=kall+allpart[i].speed[j]*allpart[i].speed[j];
        }
    }
    return kall/2.0;
}
void setv(int N,double alpha,particle* allpart){
    double normf=0;
    double normv=0;
    double temp;
    for (size_t i=0; i<N; i++) {
        for (size_t k=0; k<3; k++) {
            temp=allpart[i].force[k];
            normf=temp*temp+normf;
            temp=allpart[i].speed[k];
            normv=temp*temp+normv;
        }
    }
    normf=sqrt(normf);
    normv=sqrt(normv);
    for (size_t i=0; i<N; i++) {
        for (size_t k=0; k<3; k++) {
            allpart[i].speed[k]=(1-alpha)*allpart[i].speed[k]+alpha*normv*allpart[i].force[k]/normf;
        }
    }
}
double power(int N,particle* allpart){
    double temp=0;
    for(size_t i=0;i<N;i++){
        for(size_t k=0;k<3;k++){
            temp=temp+allpart[i].speed[k]*allpart[i].force[k];
        }
    }
    return temp;
}
void freeze(int N,particle* allpart){
    for(size_t i=0;i<N;i++){
        for(size_t k=0;k<3;k++){
            allpart[i].speed[k]=0;
        }
    }
}
void leapfrogone(int N,double Dt,double len,particle* allpart){
    double temp;
    for(size_t i=0;i<N;i++){
        for(size_t j=0;j<3;j++){
            temp=allpart[i].posit[j]+allpart[i].speed[j]*Dt+0.5*Dt*Dt*allpart[i].force[j];
            allpart[i].posit[j]=(temp/len-round(temp/len))*len;
            allpart[i].speed[j]=allpart[i].speed[j]+0.5*Dt*allpart[i].force[j];
        }
    }
}
void leapfrogtwo(int N,double Dt,particle* allpart){
    for(size_t i=0;i<N;i++){
        for(size_t j=0;j<3;j++){
            allpart[i].speed[j]=allpart[i].speed[j]+0.5*Dt*allpart[i].force[j];
        }
    }
}
//update the force for each particle.
void updateforce(int N,double len,particle* allpart){
    parnode* temp;
    double force_temp[3]={0,0,0};
    double tempdis;
    double rij;
    double dij;
    for (int i=0; i<N; i++) {
        //clear the temp force.
        for(size_t j=0;j<3;j++){
            force_temp[j]=0;
        }
        //started to sum the all force.
        temp=allpart[i].neighbor->next;
        while (temp!=NULL) {
            rij=distance(temp->index, i, len, allpart);
            dij=allpart[temp->index].radius+allpart[i].radius;
            for (size_t j=0; j<3; j++) {
                tempdis=allpart[i].posit[j]-allpart[temp->index].posit[j];
                tempdis=(tempdis/len-round(tempdis/len))*len;
                force_temp[j]=force_temp[j]+2*(1-rij/dij)*tempdis/dij/rij;
            }
            temp=temp->next;
        }
        //end sum the force.
        for(size_t j=0;j<3;j++){
            allpart[i].force[j]=force_temp[j];
        }
    }
}
//update the neighbor for one particle.
void updateallneigh(int N,int cellsize,double len,parnode* cellall[],particle* allpart){
    for (int k=0; k<N; k++) {
        updateoneneigh(k, cellsize, len, cellall, allpart);
    }
}
void updateoneneigh(int ind,int cellsize,double len,parnode* cellall[],particle* allpart){
    int box;
    double rij,dij;
    parnode* temp;
    allpart[ind].neighbor->next=NULL;
    allpart[ind].neighbor->tail=allpart[ind].neighbor;//clear the neighbor list for this particle.
    for (int i=allpart[ind].cindex[0]-2;i<allpart[ind].cindex[0]+3 ; i++) {
        for (int j=allpart[ind].cindex[1]-2; j<allpart[ind].cindex[1]+3; j++) {
            for (int k=allpart[ind].cindex[2]-2; k<allpart[ind].cindex[2]+3; k++) {
                //in these 9 cells we all potentially have iteractions with them.
                box=(i+cellsize)%cellsize+((j+cellsize)%cellsize)*cellsize+((k+cellsize)%cellsize)*cellsize*cellsize;
                temp=cellall[box]->next;
                while (temp!=NULL) {
                    if (temp->index!=ind) {
                        rij=distance(temp->index, ind, len, allpart);
                        dij=allpart[temp->index].radius+allpart[ind].radius;
                        if (rij<dij) {
                        addptolist(allpart[ind].neighbor,allpart[ind].neighbor->tail,temp->index);
                        }
                    }
                    else{
                    };
                    temp=temp->next;
                }
            }
        }
    }
};
//this function is use to update the particles in the cell.
void updatecellofp(int N,int cellsize,parnode* cellall[],particle *allpart){
    int tempindex;
    //reset the cell list
    for (size_t i=0; i<cellsize*cellsize*cellsize; i++) {
        cellall[i]->next=NULL;
        cellall[i]->tail=cellall[i];
        cellall[i]->index=-1;//use index=-1 to show that this is the head and useless.
    }
    for (int i=0; i<N; i++) {
        tempindex=allpart[i].cindex[0]+allpart[i].cindex[1]*cellsize+allpart[i].cindex[2]*cellsize*cellsize;
        addptolist(cellall[tempindex],cellall[tempindex]->tail,i);
    }
};
//this function is use to update the cell index of particles.
void updatepincell(int N,double len,double celllen,particle* allpart){
    double temp;
    for (size_t i=0; i<N;i++){
        for (size_t j=0; j<3; j++) {
            temp=allpart[i].posit[j]+0.5*len;
            allpart[i].cindex[j]=floor(temp/celllen);
        }
    }
};
void addptolist(parnode* head,parnode* tail,int in){
    parnode* temp3=(parnode*)malloc(sizeof(parnode));
    temp3->index=in;
    temp3->next=NULL;
    head->tail->next=temp3;
    head->tail=temp3;
};
//show the periodical distance of two particle.
double distance(int index1,int index2,double len,particle* allpart){
    double sum=0;
    double temp1;
    for (size_t i=0; i<3; i++) {
        temp1=allpart[index1].posit[i]-allpart[index2].posit[i];
        temp1=(temp1/len-round(temp1/len))*len;
        sum=sum+temp1*temp1;
    }
    return sqrt(sum);
};
void printplist(parnode* head){
    parnode* temp1;
    temp1=head->next;
    while (temp1!=NULL) {
        printf("%d ",temp1->index);
        temp1=temp1->next;
    }
};
