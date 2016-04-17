#ifndef TERMINAL_PRINTER_H
#define TERMINAL_PRINTER_H
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;


class TerminalPrinter{
public:
	void printStarter(){
		threadKernal();
	}
	void printStepStatus(int step,int piter ,double time,double dt,double res){
		char temp[256];
		sprintf(temp,"%15f\t%10d\t%10d\t%13.5f\t%13.5f\n",time,step,piter,dt,res);
		addMessage(string(temp));
	}
	void printSteadyStatus(int step,double res){
		char temp[256];
		sprintf(temp,"%15s\t%10d\t%10s\t%13s\t%15f\n","---",step,"---","---",res);
		addMessage(string(temp));
	}
	void printSectionHead(double timeElapse){;
		char temp[256];
		sprintf(temp,"%15s\t%10s\t%10s\t%13s\t%15s\n","TIME","CAL STEPE","ITER","DELT","MAX RES");
		addMessage(string(temp));
	}
	void printEnding(){
		char szTemp[512];
		sprintf(szTemp+strlen(szTemp), " -------------------------------------------------------------------------------- \n\n");
		sprintf(szTemp+strlen(szTemp), "                      CPU REQUIREMENTS OF NUMERICAL SOLUTION\n\n");
		sprintf(szTemp+strlen(szTemp), " -------------------------------------------------------------------------------- \n");
		sprintf(szTemp+strlen(szTemp), "    The total problem used   : %10.1f seconds of CPU time\n", 0.0);
		sprintf(szTemp+strlen(szTemp), "    includes system time     : %10.1f seconds of CPU time\n", 0.0);
		sprintf(szTemp+strlen(szTemp), "    for an average of        : %10.1f  MICROSECONDS/CELL/CYCLE.\n", 0.0);
		sprintf(szTemp+strlen(szTemp), "    Total wall clock time of : %10.1f  seconds \n", 0.0);
		addMessage(string(szTemp));
	}
	~TerminalPrinter(){
	}
	TerminalPrinter(){
	}
	void addMessage(const string& str){
		cout<<str<<endl;
	}
private:
	void threadKernal(){
		system("clear");
		printf("    1111111    111        111    1111111         111111         1111111          \n");
		printf("  11111111111  111        111  11111111111      11111111      11111111111        \n");
		printf(" 1111     1111 111        111 1111     1111    1111  1111    1111     1111       \n");
		printf(" 1111            111    111   1111            1111    1111   1111     1111       \n");
		printf(" 1111              111111     1111           1111      1111    111               \n");
		printf(" 1111               1111      1111           1111      1111      11111           \n");
		printf(" 1111               1111      1111           11111111111111         1111         \n");
		printf(" 1111               1111      1111           1111      1111  1111     1111       \n");
		printf(" 1111     1111      1111       1111    1111  1111      1111  1111     1111       \n");
		printf("  11111111111       1111        1111111111   1111      1111   11111111111        \n");
	   	printf("    1111111         1111          111111     1111      1111     1111111          \n");
		printf("\n");
		printf("\n");
    		printf("000000000000000000000000000000000000000000000000000000000000000000000000000000\n");
		printf("0 0000 0 0    00     0      0     00 0000 00     0      0 0000 0     00      0\n");
		printf("0 0000 0  000 0 00000000  000 0000 0 0000 0 00000000  000 0000 0 0000 0 000000\n");
		printf("0 0000 0 0000 00    0000  000     00 0000 0 00000000  000 0000 0     00      0\n");
		printf("0 0000 0 0000 000000 000  000 000 00 0000 0 00000000  000 0000 0 000 00 000000\n");
		printf("00    00 0000 0     0000  000 0000 00    000     000  0000    00 0000 0      0\n");
    		printf("000000000000000000000000000000000000000000000000000000000000000000000000000000\n");


		printf("                   *******************************************\n");
		printf("                     COMPUTATIONAL CODES FOR FLUID DYNAMICS\n");
		printf("                   *******************************************\n");
		printf("                   VERSION 2.0                     27 NOV 2015 \n");
		printf("                             EXECUTABLE ATTRIBUTES \n\n");
		printf("   All rights reserved. Unauthorized use, distribution or duplication is \n");
		printf("   prohibited. This product is subject to China Nuclear Power Technology \n");
		printf("   Research Institute and State Nuclear Power Software Development Center.\n");
		printf("--------------------------------------------------------------------------------\n\n");
	}
};
#endif
