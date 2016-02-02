#ifndef TERMINAL_PRINTER_H
#define TERMINAL_PRINTER_H
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <deque>
#include <queue>
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
		screenPosi.resize(11,NULL);
		vector<string> charVec;
		charVec.resize(16,NULL);
		charVec[0] = "0001111111100001110000000011100011111110000000000111111000000000011111100000000000000000000000000000";
		charVec[1] = "0011111111111001110000000011100111111111110000001111111100000011111111111000000000000000000000000000";
		charVec[2] = "0111100000111101110000000011101111000001111000011110011110000111100000111100000000000000000000000000";
		charVec[3] = "0111100000000000011100001110001111000000000000111100001111000111100000111100000000000000000000000000";
		charVec[4] = "0111100000000000000111111000001111000000000001111000000111100001110000000000000000000000000000000000";
		charVec[5] = "0111100000000000000011110000001111000000000001111000000111100000011111000000000000000000000000000000";
		charVec[6] = "0111100000000000000011110000001111000000000001111111111111100000000011110000000000000000000000000000";
		charVec[7] = "0111100000000000000011110000001111000000000001111000000111100111100000111100000000000000000000000000";
		charVec[8] = "0111100000111100000011110000000111100001111001111000000111100111100000111100000000000000000000000000";
		charVec[9] = "0011111111111000000011110000000011111111110001111000000111100011111111111000000000000000000000000000";
	   	charVec[10] = "0001111111100000000011110000000001111111000001111000000111100000011111100000000000000000000000000000";
	   //sub title in the center;
	   charVec[11] = "01000010101111001111101111110111110010000100111110111111010000101111100111111000000000000000000";
	   charVec[12] = "01000010110001010000000011000100001010000101000000001100010000101000010100000000000000000000000";
	   charVec[13] = "01000010100001001111000011000111110010000101000000001100010000101111100111111000000000000000000";
	   charVec[14] = "01000010100001000000100011000100010010000101000000001100010000101000100100000000000000000000000";
	   charVec[15] = "00111100100001011111000011000100001001111000111110001100001111001000010111111000000000000000000";
	}
	void addMessage(const string& str){
		cout<<str<<endl;
	}
private:
	vector<string> screenPosi;
	void threadKernal(){
		refressScreen();
		for(int i=0;i!=80;++i) cout<<'#';
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
	void refressScreen(){
		system("clear");
		for(vector<string>::iterator it=screenPosi.begin();it!=screenPosi.end();++it){
			int i = 0; 
			for(char* c = &it->at(0);i != 80;++i,++c){
				if(*c == '1') cout<<'#'; 
				else cout<<' ';
			}
			cout<<endl;
		}
	}
};
#endif
