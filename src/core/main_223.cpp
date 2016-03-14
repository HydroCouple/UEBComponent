//#include "uebpg.h"
//#include "nctest.h"
#include"mpi.h"
#include "uebpgdecls.h"
#include <queue>

#pragma warning(disable : 4996)
using namespace std;


int main(int argc, char* argv[])
{

	float timeControl = 0.0, timeWS = 0.0, timeSitestate = 0.0, timeTSArrays = 0.0, timeParam = 0.0, timeParamSiteInptcontrol = 0.0, timeModelRun = 0.0;
	float** outvarArray; //= new float*[70]; //[70];
	float*** aggoutvarArray;
	char conFile[256], paramFile1[256], sitevarFile1[256], inputconFile1[256], outputconFile1[256], watershedFile1[256], aggoutputconFile1[256], aggoutputFile1[256];
	char wsvarName1[256], wsycorName1[256], wsxcorName1[256];
	int **wsArray = NULL;
	int dimlen1 = 0, dimlen2 = 0, totalgrid = 0;
	int wsfillVal = -9999;
	//int  *wsmissVal = new int;
	//*wsmissVal = -9999;

	float *parvalArray = NULL;
	float SiteState[32];
	float *wsxcorArray = NULL, *wsycorArray = NULL;
	//params ParamVAlues;
	sitevar *strsvArray = new sitevar[32];
	char * svFile[32];
	char * svVarName[32];
	for (int i = 0; i < 32; i++){
		svFile[i] = new char[256];
		svVarName[i] = new char[256];
	}
	int svType[32];
	inpforcvar strinpforcArray[13];
	//outputs
	pointOutput *pOut = NULL;
	ncOutput *ncOut = NULL;
	aggOutput *aggOut = NULL;
	int npout = 0, nncout = 0, naggout = 0, nZones = 0;
	const char * zName = "Outletlocations"; //12.24.14 watershed zonning for aggregation--	
	float *z_ycor = NULL;
	float *z_xcor = NULL;
	int zoneid = 0;
	int *ZonesArr = NULL;
	float **ncoutArr;

	//inptimeseries *strintsArray[11];

	float *tcorvar[13], *tsvarArray[13], *tsvarArrayTemp[5]; // [xstride]; //assuming max nc files for a variable =5
	int ncTotaltimestep = 0;                                         //6.24.14
	int ntimesteps[5];
	int npar = 32;
	/*int NUMtimeSTEP,NREFYR,NREFMO,NREFDAY,NumOP;*/
	int ModelStartDate[3], ModelEndDate[3]; //check this 
	double ModelStartHour, ModelEndHour, ModelDt, ModelUTCOffset;
	double modelSpan;
	int numTimeStep;
	int numOut = 70;
	char headerLine[256];
	int retvalue = 0;
	int numgrid = 0;

	const char* tNameout = "time";
	int outtSteps = 0;
	int outtStride = 1, outyStep = 1, outxStep = 1;
	float* t_out;
	float out_fillVal = -9999.0;
	int outDimord = 0, aggoutDimord = 1;
	const char* tlong_name = "time";
	const char* tcalendar = "standard";
	char* uebVars[70] = { "Year", "Month", "Day", "dHour", "atff", "HRI", "Eacl", "Ema", "conZen", "Ta", "P", "V", "RH", "Qsi", "Qli", "Qnet",
		"Us", "SWE", "tausn", "Pr", "Ps", "Alb", "QHs", "QEs", "Es", "SWIT", "QMs", "Q", "FM", "Tave", "TSURFs", "cump", "cumes",
		"cumMr", "Qnet", "smelt", "refDepth", "totalRefDepth", "cf", "Taufb", "Taufd", "Qsib", "Qsid", "Taub", "Taud",
		"Qsns", "Qsnc", "Qlns", "Qlnc", "Vz", "Rkinsc", "Rkinc", "Inmax", "intc", "ieff", "Ur", "Wc", "Tc", "Tac", "QHc",
		"QEc", "Ec", "Qpc", "Qmc", "Mc", "FMc", "SWIGM", "SWISM", "SWIR", "errMB" };
	int outvarindx = 17, aggoutvarindx = 17;
	int size, rank, irank, jrank;
	double mpi_startTime = 0.0, startTimeT = 0.0, TotalTime = 0.0, totalmodelrunTime = 0.0, mpiInpreadTime = 0.0;
	double TsReadTime = 0.0, TSStartTime, ComputeStartTime, ComputeTime = 0.0, OutWriteTime;
	//mpi communication related vars
	int destRank = 0, mssgTag = 0;

	MPI::Init(argc, argv);
	//how many processes
	size = MPI::COMM_WORLD.Get_size(); //	MPI_Comm_size(MPI_COMM_WORLD,&size);
	//which rank is yours? 
	rank = MPI::COMM_WORLD.Get_rank(); //_Comm_rank(MPI_COMM_WORLD,&rank);
	//cout << "\n rank "<< rank << " of "<< size << " processes has started\n" << endl;	
	MPI_Comm myCom = MPI::COMM_WORLD;
	MPI_Info myInfo = MPI::INFO_NULL;
	if (rank == 0)
	{
		//microsecond wall time: to time block of work
		startTimeT = MPI::Wtime();
		mpi_startTime = MPI::Wtime();
		TsReadTime = 0.0;
		ComputeTime = 0.0;
	}
	//  Input Arguments		
	if (argc > 1)
	{
		//conFile = new char[sizeof(argv[0])];
		strcpy(conFile, argv[1]);
	}
	else
	{
		if (rank == 0)
			cout << "file not found exiting" << endl;
		MPI::Finalize();
		return 1;
		//cin >> conFile;
	}
	FILE* pconFile = fopen(conFile, "rt");
	fgets(headerLine, 256, pconFile);
	fscanf(pconFile, "%s\n %s\n %s\n %s\n %s\n", paramFile1, sitevarFile1, inputconFile1, outputconFile1, watershedFile1);
	fscanf(pconFile, "%s %s %s\n", wsvarName1, wsycorName1, wsxcorName1);
	fscanf(pconFile, "%s\n %s\n", aggoutputconFile1, aggoutputFile1);
	//new vs2012 appears to have issues with passing char[256] for const char*
	const char *paramFile = paramFile1, *sitevarFile = sitevarFile1, *inputconFile = inputconFile1, *outputconFile = outputconFile1, *watershedFile = watershedFile1,
		*wsvarName = wsvarName1, *wsycorName = wsycorName1, *wsxcorName = wsxcorName1, *aggoutputconFile = aggoutputconFile1, *aggoutputFile = aggoutputFile1;
	//read simulation related parameters including start and end datetimes, and model time step dt
	//#_14 tbc 5.20.13 check time format
	//fscanf(pconFile,"%f\n%f\n%f\n%f\n",startTime, endTime, dT, UTC_offset);
	fscanf(pconFile, "%d %d %d %lf\n", &ModelStartDate[0], &ModelStartDate[1], &ModelStartDate[2], &ModelStartHour);
	fscanf(pconFile, "%d %d %d %lf\n", &ModelEndDate[0], &ModelEndDate[1], &ModelEndDate[2], &ModelEndHour);
	fscanf(pconFile, "%lf\n %d %d %d\n %lf\n %d\n", &ModelDt, &outtStride, &outyStep, &outxStep, &ModelUTCOffset, &inpDailyorSubdaily);
	//close control file
	fclose(pconFile);
	//time units
	char tunits[256];
	int hhMod = (int)floor(ModelStartHour);
	int mmMod = (int)(remainder(ModelStartHour, 1.0) * 60);
	sprintf(tunits, "hours since %d-%d-%d %d:%d:00 UTC", ModelStartDate[0], ModelStartDate[1], ModelStartDate[2], hhMod, mmMod);
	const char* tUnitsout = tunits;

	//read watershed (model domain) netcdf file	
	retvalue = readwsncFile(watershedFile, wsvarName, wsycorName, wsxcorName, wsycorArray, wsxcorArray, wsArray, dimlen1, dimlen2, wsfillVal);
	//cout<<"dim1 = "<<dimlen1<<" dim2 = "<< dimlen2<<endl;
	/*printf("fillvalue= %d ",wsfillVal);
	for(int i=0;i<dimlen1;i++){
	for(int j=0;j<dimlen2;j++)
	cout<<wsArray[i][j];
	cout<<"\n";
	}*/
	//aggregation zone info
	float * wsArray1D = new float[dimlen1*dimlen2];
	for (int i = 0; i < dimlen1; i++)
	for (int j = 0; j < dimlen2; j++)
		wsArray1D[i*dimlen2 + j] = wsArray[i][j];
	std::set<int> zValues(wsArray1D, wsArray1D + (dimlen1*dimlen2));
	//cout << zValues.size() << endl;
	//std::remove_if(zValues.begin(), zValues.end(), [&wsfillVal](int a){ return a == wsfillVal; });
    std::set<int> fillSet ;//= { wsfillVal };
	//cout << "fill: " << fillSet.size() << " value: " << *(fillSet.begin())<<endl;
	std::vector<int> zVal(zValues.size());
	std::vector<int>::iterator it = std::set_difference(zValues.begin(), zValues.end(), fillSet.begin(), fillSet.end(), zVal.begin());  // exclude _FillValue
	zVal.resize(it - zVal.begin());
	//cout << zVal.size()<<endl;
	z_ycor = new float[zVal.size()];
	z_xcor = new float[zVal.size()];
	//cout << zValues.size() << endl;
	nZones = zVal.size();
	for (int iz = 0; iz < zVal.size(); iz++)
	{
		//#_12.24.14 change these with actual outlet locations coordinates
		z_ycor[iz] = 0.0;
		z_xcor[iz] = 0.0;
		//cout << zValues[iz];		
	}
	//read parameters
	//readParams(paramFile,paramValues);
	readParams(paramFile, parvalArray, npar);
	/*cout<<"param read..\n ");
	for(int i=0;i<npar;i++)
	cout<<"%f ",parvalArray[i]);	*/
	//read site vars 
	//cout<<"Reading site variable ");
	readSiteVars(sitevarFile, strsvArray); //svDefaults,svFile,svVarName,svType);
	/*cout<<"\n site variables read \n");
	for(int i=0;i<32;i++)
	cout<<"%f ",strsvArray[i].svdefValue);
	cout<<"\n");*/
	for (int i = 0; i < 32; i++)
	if (strsvArray[i].svType == 1)
	{
		//cout<<"%d %s %s\n",i, strsvArray[i].svFile,strsvArray[i].svVarName);
		retvalue = read2DNC(strsvArray[i].svFile, strsvArray[i].svVarName, strsvArray[i].svArrayValues);

		/*for(int ih=0;ih<13;ih++)
		{
		for(int jv=0;jv<16;jv++)
		cout<<"%f ",strsvArray[i].svArrayValues[ih][jv]);
		cout<<"\n");
		}*/
	}
	//read input /forcing control file--all possible entries of input control have to be provided
	readInputForcVars(inputconFile, strinpforcArray);
	modelSpan = julian(ModelEndDate[0], ModelEndDate[1], ModelEndDate[2], ModelEndHour) - julian(ModelStartDate[0], ModelStartDate[1], ModelStartDate[2], ModelStartHour);
	//model time steps
	numTimeStep = (int)ceil(modelSpan*(24 / ModelDt));
	if (rank == 0)
		cout << "number of time steps: " << " " << numTimeStep << endl;
	//allocate memory for output array	
	outvarArray = new float*[numOut];
	for (int i = 0; i < numOut; i++)
		outvarArray[i] = new float[numTimeStep];
	//total grid size to compute progress
	totalgrid = dimlen1*dimlen2;	
	//output control
	readOutputControl(outputconFile, aggoutputconFile, pOut, ncOut, aggOut, npout, nncout, naggout);
	//create output netcdf
	outtSteps = numTimeStep / outtStride;             //save SWE every outstrid'th t-step
	t_out = new float[outtSteps];
	for (int it = 0; it < outtSteps; ++it)
		t_out[it] = it*outtStride*ModelDt;      //in hours since model start time 
	ncoutArr = new float*[nncout];
	for (int icout = 0; icout < nncout; icout++)
	{
		ncoutArr[icout] = new float[outtSteps];
		for (int it = 0; it < outtSteps; ++it)
			ncoutArr[icout][it] = 0.0;
	}
	//# There should be better way than this
	aggoutvarArray = new float**[nZones];
	float * totalAgg = new float[outtSteps];
	ZonesArr = new int[nZones];
	for (int j = 0; j < nZones; j++)
	{
		ZonesArr[j] = 0;
		aggoutvarArray[j] = new float*[naggout];
		for (int i = 0; i < naggout; i++)
		{
			aggoutvarArray[j][i] = new float[outtSteps];
			for (int it = 0; it < outtSteps; it++)
				aggoutvarArray[j][i][it] = 0.0;
		}
	}
	//int numRequest = 100 * (dimlen1*dimlen2 - 1) + 13; 
	MPI::Request *isendRequest = new MPI::Request[nncout];     // , ireceiveRquest, frontRequest, backRquest;
	MPI::Status isendStatus;
	//std::queue<MPI::Request> requestQueue;
	//pushback -add to que;    front --get first; pop --revmove from queue .empty check queue empty
	//wait for all processes to complete
	//MPI::COMM_WORLD.Barrier();
	if (rank == 0)
	{
		mpiInpreadTime = MPI::Wtime() - mpi_startTime;
		mpi_startTime = MPI::Wtime();
		for (int icout = 0; icout < nncout; icout++)
			retvalue = create3DNC_uebOutputs(ncOut[icout].outfName, (const char*)ncOut[icout].symbol, (const char*)ncOut[icout].units, tNameout, tUnitsout,
			tlong_name, tcalendar, outtSteps, outDimord, t_out, &out_fillVal, watershedFile, wsvarName, wsycorName, wsxcorName);
		//create aggregate ouput file
		retvalue = create3DNC_uebAggregatedOutputs(aggoutputFile, aggOut, naggout, tNameout, tUnitsout, tlong_name, tcalendar, outtSteps, aggoutDimord, t_out, &out_fillVal,
			watershedFile, wsvarName, wsycorName, wsxcorName, nZones, zName, z_ycor, z_xcor);
		destRank = 0;		
		//**** Space loop starts here	
		for (irank = 0; irank < dimlen1; irank++)
		{
			for (jrank = 0; jrank < dimlen2; jrank++)
			{
				destRank++;           //rank of destination process
				if (destRank >(size - 1))
					destRank = 1;
				//if ((240 < irank && irank < 244) && (222 < jrank && jrank < 226))		                              //( fmod(float(irank),50)==0  && fmod(float(jrank),60) ==0
				{
					//cout <<"cell val: "<< wsArray[irank][jrank] << endl;
					//cout << "fill val: "<<wsfillVal << endl;
					if (wsArray[irank][jrank] != wsfillVal)      //Fill value in ws => non-computation cell, exit loop
					{
						//read time series forcing data
						if (strsvArray[16].svType != 3)  //no accumulation zone 
							//***tbc what happens if it is accumulation zone?
						{
							for (int it = 0; it < 13; it++)            //it < 13,          12.18.14
							{
								if (strinpforcArray[it].infType == 0)
									readTextData(strinpforcArray[it].infFile, tsvarArray[it], ncTotaltimestep);
								else if (strinpforcArray[it].infType == 1)     // == 0
								{
									ncTotaltimestep = 0;
									for (int numNc = 0; numNc < strinpforcArray[it].numNcfiles; numNc++) //if multiple netcdf for a single variable, they are read one by one and copied to single array
									{
										//read 3D netcdf (regridded array processed by uebInputs)
										char numtoStr[256];
										sprintf(numtoStr, "%d", numNc);
										char tsInputfile[256];
										strcpy(tsInputfile, strinpforcArray[it].infFile);
										strcat(tsInputfile, numtoStr);
										strcat(tsInputfile, ".nc");
										//cout<<"%s\n",tsInputfile);
										retvalue = readNC_TS(tsInputfile, strinpforcArray[it].infvarName, strinpforcArray[it].inftimeVar,
											wsycorName, wsxcorName, tsvarArrayTemp[numNc], tcorvar[it], irank, jrank, ntimesteps[numNc]);
										ncTotaltimestep += ntimesteps[numNc];
										/*for(int tps=0;tps<ncTotaltimestep;tps++)
											cout << "  " << tsvarArrayTemp[numNc][xstrt][tps];
											cout<<"  "<<ncTotaltimestep<<endl;*/
									}
									tsvarArray[it] = new float[ncTotaltimestep];
									int tinitTime = 0;
									for (int numNc = 0; numNc < strinpforcArray[it].numNcfiles; numNc++)
									{
										for (int tts = 0; tts < ntimesteps[numNc]; tts++)
											tsvarArray[it][tts + tinitTime] = tsvarArrayTemp[numNc][tts];
										tinitTime += ntimesteps[numNc];
									}
									/*for (int tps = 0; tps < tinitTime; tps++)
									cout << "  " << tsvarArray[xstrt][it][tps] << " ";
									cout << endl << " " << tinitTime << endl;*/
								}
								else       //if (strinpforcArray[it].infType == 2 || strinpforcArray[it].infType == -1)
								{
									//######TBC 6.20.13 better way to handle this is needed
									tsvarArray[it] = new float[2];
									ncTotaltimestep = 2;
									//just copy the default value if a single value is the option									
									tsvarArray[it][0] = strinpforcArray[it].infType;
									tsvarArray[it][1] = strinpforcArray[it].infdefValue;
								}
								//send data
								//cmmented out 2.4.15 destRank = (irank * (size-1)) / dimlen1 + 1;   //distribute evenly then add the extra cells at the beginning 2.3.14
								if (size > 1)
								{
									mssgTag = 100 * (irank*dimlen2 + jrank) + it;            //message tag 
									MPI::COMM_WORLD.Send(tsvarArray[it], ncTotaltimestep, MPI::FLOAT, destRank, mssgTag);
								}
							}
							/*int stepcontrol = jrank + xstrt;
							cout << "jrank = " << stepcontrol << " ";*/
							/*cout << endl << "Tmin " << endl;
							for (int tps = 0; tps < 100; tps++)
							cout << "  " << tsvarArray[xstrt][10][tps] << " ";
							cout << endl << "Tmax " << endl;
							for (int tps = 0; tps < 100; tps++)
							cout << "  " << tsvarArray[xstrt][11][tps] << " ";
							cout << endl << "Prec " << endl;
							for (int tps = 0; tps < 100; tps++)
							cout << "  " << tsvarArray[xstrt][1][tps] << " ";
							cout << endl << "Wind " << endl;
							for (int tps = 0; tps < 100; tps++)
							cout << "  " << tsvarArray[xstrt][2][tps] << " ";*/
							//write nc outputs
							if (size > 1)
							{
								/*for (int icout = 0; icout < nncout; icout++)
								{
									mssgTag = 100 * (irank*dimlen2 + jrank) + icout;            //message tag 
									backRquest = MPI::COMM_WORLD.Irecv(t_out, outtSteps, MPI::FLOAT, destRank, mssgTag);
									requestQueue.push(backRquest);									
								}*/
							}
							else //size = 1
							{
								//call runUEB ---	
								RUNUEB(tsvarArray, SiteState, parvalArray, outvarArray, ModelStartDate, ModelStartHour, ModelEndDate, ModelEndHour, ModelDt, ModelUTCOffset);
								//write nc outputs
								for (int icout = 0; icout < nncout; icout++)
								{
									for (int vindx = 0; vindx < 70; vindx++)
									{
										if (strcmp(ncOut[icout].symbol, uebVars[vindx]) == 0)
										{
											outvarindx = vindx;
											break;
										}
									}
									for (int it = 0; it < outtSteps; ++it)
										t_out[it] = outvarArray[outvarindx][outtStride*it];         //use timeStiride to sample outputs if it is dense (e.g hourly data for a year may be too big to save in one nc file)
									//write var values
									retvalue = WriteTSto3DNC((const char*)ncOut[icout].outfName, (const char*)ncOut[icout].symbol, outDimord, irank, jrank, outtSteps, t_out);
								}
								//point outputs
								for (int ipout = 0; ipout < npout; ipout++)
								{
									if (irank == pOut[ipout].ycoord && (jrank) == pOut[ipout].xcoord)
									{
										FILE* pointoutFile = fopen((const char*)pOut[ipout].outfName, "w");
										for (int istep = 0; istep < numTimeStep; istep++)
										{
											fprintf(pointoutFile, "\n %d %d %d %8.3f ", (int)outvarArray[0][istep], (int)outvarArray[1][istep], (int)outvarArray[2][istep], outvarArray[3][istep]);
											for (int vnum = 4; vnum < 70; vnum++)
												fprintf(pointoutFile, " %16.6f ", outvarArray[vnum][istep]);
										}
										fclose(pointoutFile);
									}
								}
								//#_??aggregated outputs 12.24.14
								zoneid = wsArray[irank][jrank] - 1;
								ZonesArr[zoneid] += 1;
								for (int iagout = 0; iagout < naggout; iagout++)
								{
									for (int vindx = 0; vindx < 70; vindx++)
									{
										if (strcmp(aggOut[iagout].symbol, uebVars[vindx]) == 0)
										{
											aggoutvarindx = vindx;
											break;
										}
									}
									for (int it = 0; it < outtSteps; it++)
										aggoutvarArray[zoneid][iagout][it] += outvarArray[aggoutvarindx][outtStride*it];
								}
							}
						}//end if(! accumulation zone)
					}  //if if (wsArray[irank][jrank] !=  (int)wsFillval)			
				}  //if( fmod(irank,50)==0  && fmod(jrank,60) ==0)
			}	//for(jrank=0; jrank < dimlen2; jrank++)	
			//grid count progress is calculated and written here
			if (size == 1)
			{
				numgrid += dimlen2;   // *size;
				cout << "\r   percent completed: " << ((float)numgrid / (float)totalgrid)*100.0 << " %" << endl;
				fflush(stdout);
			}
		} //end dim1
		destRank = 0;
		if (size > 1)
		{
			for (irank = 0; irank < dimlen1; irank++)
			{
				for (jrank = 0; jrank < dimlen2; jrank++)
				{
					destRank++;           //rank of destination process
					if (destRank >(size - 1))
						destRank = 1;
					//if ((240 < irank && irank < 244) && (222 < jrank && jrank < 226))		                              //( fmod(float(irank),50)==0  && fmod(float(jrank),60) ==0
					{
						if (wsArray[irank][jrank] != wsfillVal)      //Fill value in ws => non-computation cell, exit loop
						{
							if (strsvArray[16].svType != 3)  //no accumulation zone 							
							{
								for (int icout = 0; icout < nncout; icout++)
								{
									mssgTag = 100 * (irank*dimlen2 + jrank) + icout;            //message tag									
									MPI::COMM_WORLD.Recv(t_out, outtSteps, MPI::FLOAT, destRank, mssgTag);
									//write var values
									retvalue = WriteTSto3DNC((const char*)ncOut[icout].outfName, (const char*)ncOut[icout].symbol, outDimord, irank, jrank, outtSteps, t_out);
								}
							}//end if(! accumulation zone)
						}  //if if (wsArray[irank][jrank] !=  (int)wsFillval)			
					}  //if( fmod(irank,50)==0  && fmod(jrank,60) ==0)
				}	//for(jrank=0; jrank < dimlen2; jrank++)
				//grid count progress is calculated and written here						
				numgrid += dimlen2;        // *size;				
				cout << "\r  Writing output, percent completed: " << ((float)numgrid / (float)totalgrid)*100.0 << " %" << endl;
				fflush(stdout);
			} //end dim1
		}
	}
	else //rank !=0
	{
		//receive message			
		for (irank = 0; irank < dimlen1; irank++)
		{
			for (jrank = (rank - 1); jrank < dimlen2; jrank += (size - 1))
			{
				//if ((240 < irank && irank < 244) && (222 < jrank && jrank < 226))		                     
				{
					//cout <<"cell val: "<< wsArray[irank][jrank] << endl;
					//cout << "fill val: "<<wsfillVal << endl;
					if (wsArray[irank][jrank] != wsfillVal)      //Fill value in ws => non-computation cell, exit loop
					{
						//track grid cell
						uebCellY = irank;
						uebCellX = jrank;
						for (int is = 0; is < 32; is++)
						{
							if (strsvArray[is].svType == 1)
								//#_TBC 6.4.13 cast may not be necessariy 
								SiteState[is] = strsvArray[is].svArrayValues[irank][jrank];
							else
								SiteState[is] = strsvArray[is].svdefValue;
						}
						//read time series forcing data
						if (strsvArray[16].svType != 3)  //no accumulation zone 
							//***tbc what happens if it is accumulation zone?
						{
							for (int it = 0; it < 13; it++)            //it < 13,          12.18.14
							{
								//receive data 
								mssgTag = 100 * (irank*dimlen2 + jrank) + it;            //message tag from 0 to dimlen1*dimlen2-1
								tsvarArray[it] = new float[numTimeStep + 100];                 //+100 add extra space to be safe 2.4.15
								MPI::COMM_WORLD.Recv(tsvarArray[it], numTimeStep + 100, MPI::FLOAT, 0, mssgTag);
							}
							//call runUEB ---	
							RUNUEB(tsvarArray, SiteState, parvalArray, outvarArray, ModelStartDate, ModelStartHour, ModelEndDate, ModelEndHour, ModelDt, ModelUTCOffset);
							//write nc outputs
							for (int icout = 0; icout < nncout; icout++)
							{
								for (int vindx = 0; vindx < 70; vindx++)
								{
									if (strcmp(ncOut[icout].symbol, uebVars[vindx]) == 0)
									{
										outvarindx = vindx;
										break;
									}
								}
								for (int it = 0; it < outtSteps; ++it)
									ncoutArr[icout][it] = outvarArray[outvarindx][outtStride*it];         //use timeStiride to sample outputs if it is dense (e.g hourly data for a year may be too big to save in one nc file)								
								mssgTag = 100*(irank*dimlen2 + jrank) + icout;            //message tag 
								//initiate output send 2.5.15
								isendRequest[icout] =  MPI::COMM_WORLD.Isend(ncoutArr[icout],outtSteps, MPI::FLOAT, 0, mssgTag);								
							}
							//point outputs
							for (int ipout = 0; ipout < npout; ipout++)
							{
								if (irank == pOut[ipout].ycoord && (jrank) == pOut[ipout].xcoord)
								{
									FILE* pointoutFile = fopen((const char*)pOut[ipout].outfName, "w");
									for (int istep = 0; istep < numTimeStep; istep++)
									{
										fprintf(pointoutFile, "\n %d %d %d %8.3f ", (int)outvarArray[0][istep], (int)outvarArray[1][istep], (int)outvarArray[2][istep], outvarArray[3][istep]);
										for (int vnum = 4; vnum < 70; vnum++)
											fprintf(pointoutFile, " %16.6f ", outvarArray[vnum][istep]);
									}
									fclose(pointoutFile);
								}
							}
							//#_??aggregated outputs 12.24.14
							zoneid = wsArray[irank][jrank] - 1;
							ZonesArr[zoneid] += 1;
							for (int iagout = 0; iagout < naggout; iagout++)
							{
								for (int vindx = 0; vindx < 70; vindx++)
								{
									if (strcmp(aggOut[iagout].symbol, uebVars[vindx]) == 0)
									{
										aggoutvarindx = vindx;
										break;
									}
								}
								for (int it = 0; it < outtSteps; it++)
									aggoutvarArray[zoneid][iagout][it] += outvarArray[aggoutvarindx][outtStride*it];
							}
							//finish/wait for send to finsih
							for (int icout = 0; icout < nncout; icout++)
							{
								isendRequest[icout].Wait(isendStatus);
								//if (isendStatus.Get_error() !=)       //2.5.15 to do check sedn status fo rerror
							}
							//debug outputs
							/*if (irank % outyStep == 0 && (jrank + xstrt) % outxStep == 0)
							{
							char testPrint[256];
							char ind[256];
							strcpy(testPrint, "ZTest");
							sprintf(ind, "%d", irank);
							strcat(testPrint, ind);
							strcat(testPrint, "_");
							sprintf(ind, "%d", jrank + xstrt);
							strcat(testPrint, ind);
							strcat(testPrint, ".txt");
							FILE* testoutFile = fopen(testPrint, "w");
							for (int istep = 0; istep < numTimeStep; istep++)
							{
							fprintf(testoutFile, "\n %d %d %d %8.3f ", (int)outvarArray[0][istep], (int)outvarArray[1][istep], (int)outvarArray[2][istep], outvarArray[3][istep]);
							for (int vnum = 4; vnum < 70; vnum++)
							fprintf(testoutFile, " %16.6f ", outvarArray[vnum][istep]);
							}
							fclose(testoutFile);
							}*/
						}//end if(! accumulation zone)
					}  //if if (wsArray[irank][jrank] !=  (int)wsFillval)			
				}  //if( fmod(irank,50)==0  && fmod(jrank,60) ==0)
			}	//for(jrank=0; jrank < dimlen2; jrank++)		
			//grid count progress is calculated and written here
			if (rank == 1)
			{
				numgrid += dimlen2;        // *size;				
				cout << "\r  Computation, percent completed: " << ((float)numgrid / (float)totalgrid)*100.0 << " %" << endl;
				fflush(stdout);
			}
		} //end dim1
	}   //if rank !=0
	MPI::COMM_WORLD.Barrier();
	//aggregation/ reduction 
	int rankrec = 0;               //receiver rank 
	int totalZonecells = 1;
	for (int izone = 0; izone < nZones; izone++)
	{
		rankrec = izone*size / nZones;
		for (int iagout = 0; iagout < naggout; iagout++)
		{
			MPI::COMM_WORLD.Reduce(aggoutvarArray[izone][iagout], totalAgg, outtSteps, MPI::FLOAT, MPI::SUM, rankrec);
			MPI::COMM_WORLD.Reduce(&ZonesArr[izone], &totalZonecells, 1, MPI::INT, MPI::SUM, rankrec);
			//#_12.28.14 aggregation operation needs defining
			if (rank == rankrec)
			{
				if (strcmp(aggOut[iagout].aggop, "AVE") == 0)
				for (int it = 0; it < outtSteps; it++)
					totalAgg[it] /= totalZonecells;
				retvalue = Write_uebaggTS_toNC(aggoutputFile, aggOut[iagout].symbol, aggoutDimord, izone, outtSteps, totalAgg);
			}
		}
	}
	/*cout<<"Process "<<rank<<" finished computation"<<endl;
	fflush(stdout);*/
	//deallocate memory ====#_*_#______Needs revisiting; some of the arrays are not deleted 6.23.13
	for (int i = 0; i < dimlen1; i++)
		delete[] wsArray[i];
	delete[] wsArray;
	//delete[] wsycorArray;
	//delete[] wsxcorArray;
	delete[] parvalArray;
	for (int i = 0; i < 32; i++)
	{
		if (strsvArray[i].svType == 1)
		{
			for (int j = 0; j < dimlen1; j++)
				delete[] strsvArray[i].svArrayValues[j];
			delete[] strsvArray[i].svArrayValues;
		}
	}
	delete[] strsvArray;
	//delete[] tsvarArray[kx];
	for (int it = 0; it < 13; it++)       //10-->12   6.26.14
	{
		delete[] tsvarArray[it];
	}
	//delete[] tsvarArray;
	/*for (int it = 0; it < 5; it++)       //10-->12   6.26.14
	{
		delete[] tsvarArrayTemp[it];
	}*/
	//delete[] tsvarArrayTemp;

	//delete[] tcorvar;
	for (int it = 0; it < 13; it++)
	{
		if (strinpforcArray[it].infType == 1)
			delete[] tcorvar[it];
	}
	for (int zk = 0; zk < nZones; zk++)
	{
		for (int ig = 0; ig < naggout; ig++)
			delete[] aggoutvarArray[zk][ig];
		delete[] aggoutvarArray[zk];
	}
	delete[] aggoutvarArray;
	/*for(int k=0 ;k<numOut; k++)
		delete[] outvarArray[k];
		delete []outvarArray; 	*/
	cout << "Process " << rank << " finished" << endl;
	fflush(stdout);
	//wait for all processes to finish
	//MPI::COMM_WORLD.Barrier();

	if (rank == 0)
	{
		TotalTime = MPI::Wtime() - startTimeT;                   //(float)1000*(endTimeT - startTimeT)/CLOCKS_PER_SEC;
		cout << "\n Time in micro seconds" << endl;
		cout << "Reading param  site state input control:  " << mpiInpreadTime << endl;
		cout << "Reading input TS arrays:  " << TsReadTime << endl;
		cout << "Model simulation run time:  " << ComputeTime << endl;
		//cout<<"Write outputs time: "<<OutWriteTime<<endl;
		cout << "Total time of including overhead :  " << TotalTime << endl;
		cout << "Done! return value: " << retvalue << endl;
		fflush(stdout);
	}

exitlab:
	MPI::Finalize();
	getchar();
	return 0;
}


