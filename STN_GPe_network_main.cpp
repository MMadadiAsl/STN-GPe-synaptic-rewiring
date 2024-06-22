
//********************************************************************************************
//********************************************************************************************

	#include <iostream>
	#include <cmath>
	#include <fstream>
	#include <vector>
	#include "random_number.cpp"
	#include "spike_statistics.cpp"

	using namespace std;

//********************************************************************************************
//********************************************************************************************

int main(){

	// Define the size of each population
	const int GPe_Num = 460;
	const int STN_Num = 136;

	const float PI = 4.0 * atan(1.0);	// Definition of the pi number.
	long iseed = -80L;			// Seed for random numbers.

	float P_GG = 0.05;			// GPe-GPe connection probability.
	float P_GS = 0.02;			// GPe-STN connection probability.
	float P_SG = 0.05;			// STN-GPe connection probability.

	float dt = 0.02;			// Time step of the simulation.
	int T = 1250000;			// Total steps of the simulation.

	// Define the variables
	float **VGPe, **VSTN, **sS, **sG, *Gmean, **nG, **hG, **rG, **CaG, **nS, **hS, **rS, **CaS;
	float ILG, IKG, INaG, ITG, ICaG, IAHPG, ILS, IKS, INaS, ITS, ICaS, IAHPS;
	float ainfG[GPe_Num], sinfG[GPe_Num], rinfG[GPe_Num], minfG[GPe_Num], ninfG[GPe_Num], hinfG[GPe_Num], HinfG[GPe_Num], taunG[GPe_Num], tauhG[GPe_Num];
	float ainfS[STN_Num], sinfS[STN_Num], rinfS[STN_Num], minfS[STN_Num], ninfS[STN_Num], hinfS[STN_Num], HinfS[STN_Num], taunS[STN_Num], tauhS[STN_Num];
	float binfS[STN_Num], taurS[STN_Num], IappG[GPe_Num], IappS[STN_Num];
	int tf_GPe[GPe_Num], tf_STN[STN_Num];

	VGPe = new float* [GPe_Num];
	sG = new float* [GPe_Num];
	nG = new float* [GPe_Num];
	hG = new float* [GPe_Num];
	rG = new float* [GPe_Num];
	CaG = new float* [GPe_Num];

	VSTN = new float* [STN_Num];
	sS = new float* [STN_Num];
	nS = new float* [STN_Num];
	hS = new float* [STN_Num];
	rS = new float* [STN_Num];
	CaS = new float* [STN_Num];
	Gmean = new float [T];

	for (int i = 0; i < GPe_Num; i++){
		VGPe[i] = new float [T];
		sG[i] = new float [T];
		nG[i] = new float [T];
		hG[i] = new float [T];
		rG[i] = new float [T];
		CaG[i] = new float [T];
	}

	for (int i = 0; i < STN_Num; i++){
		VSTN[i] = new float [T];
		sS[i] = new float [T];
		nS[i] = new float [T];
		hS[i] = new float [T];
		rS[i] = new float [T];
		CaS[i] = new float [T];
	}

//********************************************************************************************
//********************************************************************************************

	// GPe parameters
	float VNaG = 55.0, VKG = -80.0, VCaG = 120.0, VLG = -55.0, VsynSG = 0.0, VsynSS = 15.0;
	float gNaG = 120.0, gKG = 30.0, gAHPG = 30.0, gTG = 0.5, gCaG = 0.15, gLG = 0.1;
	float thaG = -57.0, sigaG = 2.0, thsG = -35.0, sigsG = 2.0, thrG = -70.0, sigrG = -2.0, taurG = 30.0;
	float thmG = -37.0, sigmG = 10.0, thnG = -50.0, signG = 14.0;
	float taun0G = 0.05, taun1G = 0.27, thnGt = -40.0, snG = -12.0, thHG = -57.0, sigHG = 2.0;;
	float thhG = -58.0, sighG = -12.0, tauh0G = 0.05, tauh1G = 0.27, thhGt = -40.0, shG = -12.0;
	float k1G = 30.0, kCaG = 20.0, epsG = 0.0001;
	float phirG = 1.0, phinG = 0.05, phihG = 0.05, Cm = 1.0;
	float alphaG = 2.0, betaG = 0.08, thG = 20.0;

	// STN parameters
	float VNaS = 55.0, VKS = -80.0, VCaS = 140.0, VLS = -60.0, VsynGS = -85.0, VsynGG = -100.0;
	float gNaS = 37.5, gKS = 45.0, gAHPS = 9.0, gTS = 0.5, gCaS = 0.5, gLS = 2.25;
	float thaS = -63.0, sigaS = 7.8, thsS = -39.0, sigsS = 8.0, thrS = -67.0, sigrS = -2.0;
	float thmS = -30.0, sigmS = 15.0, thnS = -32.0, signS = 8.0, thetab = 0.4, sigmab = -0.1;
	float taun0S = 1.0, taun1S = 100.0, thnSt = -80.0, snS = -26.0, thHS = -39.0, sigHS = 8.0;
	float thhS = -39.0, sighS = -3.1, tauh0S = 1.0, tauh1S = 500, thhSt = -57.0, shS = -3.0;
	float srS = -2.2, thrSt = 68.0, taur0S = 40.0, taur1S = 17.5;
	float k1S = 15.0, kCaS = 22.5, epsS = 0.00005;
	float phirS = 0.2, phinS = 0.75, phihS = 0.75;
	float alphaS = 5.0, betaS = 1.0, thS = 30.0; 

	// STDP parameters
	float eta_STDP = 0.005, tau_STDP = 20.0, alpha_STDP = 0.2;

	// Other parameters
	float gmax = 0.5, gmin = 0.05, delta_g = 0.0, LFP_GPe = 0.0, LFP_STN = 0.0, mean_sum = 0.0, mean_rate_GPe, mean_rate_STN;
	float IsynGS = 0.0, IsynSG = 0.0, IsynGG = 0.0, IsynG = 0.0, IsynS = 0.0, random_num, delta_t, rate_GPe[GPe_Num], rate_STN[STN_Num];
	float X_GPe, X_STN, X_syn, variance_GPe[GPe_Num], variance_STN[STN_Num], VGPe_dummy[int(0.1 * T)], VSTN_dummy[int(0.1 * T)], tp_GPe[5], tp_STN[5];
	float variance_global_GPe, variance_global_STN, Global_V_GPe[int(0.1 * T)], Global_V_STN[int(0.1 * T)];
	float amp_stim = 100.0, delta_stim = 5.0, freq_stim = 15.0, time_shift = 30.0 * dt, I_stim_GPe, I_stim_STN, tn_GPe, tn_STN, T_OFF = 500.0;
	int counterGG = 0, counterSG = 0, counterGS = 0, GPe_spike_counter = 0, STN_spike_counter = 0, counter = 0, counter3 = 0, counter_mean = 0, n_GPe = 0, n_STN = 0;
	int num_spike_GPe[GPe_Num] = {0}, num_spike_STN[STN_Num] = {0}, counter4 = 0, counter5 = 0, nb_GPe = 0, np_GPe = 0, nb_STN = 0, np_STN = 0, p = 5;

	for (int i = 0; i < GPe_Num; i++){
		VGPe[i][0] = 1.0;
		sG[i][0] = 0.1;
		nG[i][0] = 0.1;
		hG[i][0] = 0.8;
		rG[i][0] = 0.1;
		CaG[i][0] = 0.1;
		tf_GPe[i] = 0;
		IappG[i] = (0.02 * gasdev(&iseed)) + (-2.9);
		//IappG[i] = (0.02 * gasdev(&iseed)) + (-0.1);
	}

	for (int i = 0; i < STN_Num; i++){
		VSTN[i][0] = 1.0;
		sS[i][0] = 0.1;
		nS[i][0] = 0.1;
		hS[i][0] = 0.1;
		rS[i][0] = 0.1;
		CaS[i][0] = 0.1;
		tf_STN[i] = 0;
		IappS[i] = (0.02 * gasdev(&iseed)) + (6.9);
		//IappS[i] = (0.02 * gasdev(&iseed)) + (0.5);
	}

	// Write to a text file
	ofstream outfile1 ("structural_connectivity.txt");
	ofstream outfile2 ("mean_coupling.txt");
	ofstream outfile3 ("raster_GPe.txt");
	ofstream outfile4 ("raster_STN.txt");
	ofstream outfile5 ("LFP_GPe.txt");
	ofstream outfile6 ("LFP_STN.txt");

//********************************************************************************************
//********************************************************************************************

	// Initializing structural connectivity network
	// Create an adjacency and strength matrix
	float adj_matrix[STN_Num + GPe_Num][STN_Num + GPe_Num], strength_matrix[STN_Num + GPe_Num][STN_Num + GPe_Num];

	// Populate the adjacency matrix with intra-population connectivity: STN --> STN
        for (int i = 0; i < STN_Num; i++){
		for (int j = 0; j < STN_Num; j++){
			adj_matrix[i][j] = 0;
			strength_matrix[i][j] = 0;
		}
	}

	// Populate the adjacency matrix with intra-population connectivity: GPe --> GPe
	for (int i = STN_Num; i < STN_Num + GPe_Num; i++){
		for (int j = STN_Num; j < STN_Num + GPe_Num; j++){
			if (i == j){
				adj_matrix[i][j] = 0;
				strength_matrix[i][j] = 0;
			}
			else{
				random_num = ran2(&iseed);
				if (random_num <= P_GG){
					adj_matrix[i][j] = -1;
					strength_matrix[i][j] = (gasdev(&iseed) * 0.03) + (gmax / 2.0);
				}
				else{
					adj_matrix[i][j] = 0;
					strength_matrix[i][j] = 0;
				}
			}
		}
	}

	// Populate the adjacency matrix with inter-population connectivity: GPe --> STN
	for (int i = STN_Num; i < STN_Num + GPe_Num; i++){
		for (int j = 0; j < STN_Num; j++){
			random_num = ran2(&iseed);
			if (random_num <= P_GS){
				adj_matrix[i][j] = -1;
				strength_matrix[i][j] = (gasdev(&iseed) * 0.03) + (gmax / 2.0);
				//if (strength_matrix[i][j] <= 0.10) adj_matrix[i][j] = 0;
				//else adj_matrix[i][j] = -1;
			}
			else{
				adj_matrix[i][j] = 0;
				strength_matrix[i][j] = 0;
			}
		}
	}

	// Populate the adjacency matrix with inter-population connectivity: STN --> GPe
	for (int i = 0; i < STN_Num; i++){
		for (int j = STN_Num; j < STN_Num + GPe_Num; j++){
			random_num = ran2(&iseed);
			if (random_num <= P_SG){
				adj_matrix[i][j] = 1;
				strength_matrix[i][j] = (gasdev(&iseed) * 0.03) + (gmax / 2.0);
			}
			else{
				adj_matrix[i][j] = 0;
				strength_matrix[i][j] = 0;
			}
		}
	}

	// Print the adjacency matrix
	//for (int i = 0; i < STN_Num + GPe_Num; i++) for (int j = 0; j < STN_Num + GPe_Num; j++) outfile1 << i << '\t' << j << '\t' << adj_matrix[i][j] << endl;
	//for (int i = 0; i < STN_Num + GPe_Num; i++) for (int j = 0; j < STN_Num + GPe_Num; j++) outfile1 << i << '\t' << j << '\t' << strength_matrix[i][j] << endl;

//**************************************************************************************************
//**************************************************************************************************

	// Time loop
	for (int t = 0; t < T - 1; t++){

		//Mean coupling
		for (int i = STN_Num; i < STN_Num + GPe_Num; i++){ 
			for (int j = 0; j < STN_Num; j++){
				if (adj_matrix[i][j] == -1){
				mean_sum = mean_sum + strength_matrix[i][j];
				counter_mean++;
				}
			}
		}

		Gmean[t] = (1. / float(counter_mean)) * mean_sum; 
		outfile2 << t * dt / 1000. << '\t' << Gmean[t] << endl;
		mean_sum = 0.0;
		counter_mean = 0;

//**************************************************************************************************

			if (t > 250000 && t < 500000){

				// Stimulation current: Monopolar - continuous
				//I_stim_GPe = amp_stim * Step_Function(sin(2.0 * PI * t * dt * (freq_stim * 0.001))) * (1.0 - Step_Function(sin(2.0 * PI * (t + delta_stim) * dt * (freq_stim * 0.001))));
				
				// Stimulation current: Monopolar - burst
				tp_GPe[np_GPe] = ((1000.0 * np_GPe * (1.0 / freq_stim)) + (nb_GPe * (T_OFF + 291.7))) / dt;
				if (tp_GPe[np_GPe] <= t && t <= (tp_GPe[np_GPe] + delta_stim)) I_stim_GPe = amp_stim;
				else I_stim_GPe = 0.0;
				if (t > (tp_GPe[np_GPe] + delta_stim)) np_GPe++;
				if (np_GPe > p){
					np_GPe = 0;
					nb_GPe++;
				}

				// Stimulation current: Biopolar - continuous
				//tn_GPe = ((1000.0 * n_GPe) / freq_stim) / dt;
				//if (tn_GPe <= t && t < (tn_GPe + delta_stim)) I_stim_GPe = amp_stim;
				//else if ((tn_GPe + delta_stim) <= t && t < (tn_GPe + (2.0 * delta_stim))) I_stim_GPe = - amp_stim;
				//else I_stim_GPe = 0.0; 
				//if (t > tn_GPe + (2.0 * delta_stim)) n_GPe++;
				
			}

			else I_stim_GPe = 0.0;
			//outfile7 << t * dt / 1000.0 << '\t' << I_stim_GPe << endl;

			if (t > (250000 + time_shift) && t < (500000 + time_shift)){
				
				// Stimulation current: Monopolar - continuous
				//I_stim_STN = amp_stim * Step_Function(sin(2.0 * PI * t * dt * (freq_stim * 0.001))) * (1.0 - Step_Function(sin(2.0 * PI * (t + delta_stim) * dt * (freq_stim * 0.001))));
				
				// Stimulation current: Monopolar - burst
				tp_STN[np_STN] = ((1000.0 * np_STN * (1.0 / freq_stim)) + (nb_STN * (T_OFF + 291.7))) / dt;
				if (tp_STN[np_STN] <= t && t <= (tp_STN[np_STN] + delta_stim)) I_stim_STN = amp_stim;
				else I_stim_STN = 0.0;
				if (t > (tp_STN[np_STN] + delta_stim)) np_STN++;
				if (np_STN > p){
					np_STN = 0;
					nb_STN++;
				}
				
				// Stimulation current: Biopolar - continuous
				//tn_STN = ((1000.0 * n_STN) / freq_stim) / dt;
				//if (tn_STN <= t && t < (tn_STN + delta_stim)) I_stim_STN = amp_stim;
				//else if ((tn_STN + delta_stim) <= t && t < (tn_STN + (2.0 * delta_stim))) I_stim_STN = - amp_stim;
				//else I_stim_STN = 0.0; 
				//if (t > tn_STN + (2.0 * delta_stim)) n_STN++;
				
			}

			else I_stim_STN = 0.0;
			//outfile8 << t * dt / 1000.0 << '\t' << I_stim_STN << endl;


		// GPe Cell
		for (int i = 0; i < GPe_Num; i++){
			ainfG[i] = 1.0 / (1.0 + exp(-(VGPe[i][t] - thaG) / sigaG));
			sinfG[i] = 1.0 / (1.0 + exp(-(VGPe[i][t] - thsG) / sigsG));
			rinfG[i] = 1.0 / (1.0 + exp(-(VGPe[i][t] - thrG) / sigrG));
			minfG[i] = 1.0 / (1.0 + exp(-(VGPe[i][t] - thmG) / sigmG));
			ninfG[i] = 1.0 / (1.0 + exp(-(VGPe[i][t] - thnG) / signG));
			hinfG[i] = 1.0 / (1.0 + exp(-(VGPe[i][t] - thhG) / sighG));
			HinfG[i] = 1.0 / (1.0 + exp(-(VGPe[i][t] - thG - thHG) / sigHG));
			taunG[i] = taun0G + taun1G / (1.0 + exp(-(VGPe[i][t] - thnGt) / snG));
			tauhG[i] = tauh0G + tauh1G / (1.0 + exp(-(VGPe[i][t] - thhGt) / shG));
	
			ILG = gLG * (VGPe[i][t] - VLG);
			IKG = gKG * nG[i][t]*nG[i][t]*nG[i][t]*nG[i][t] * (VGPe[i][t] - VKG);
			INaG = gNaG * minfG[i]*minfG[i]*minfG[i]*hG[i][t] * (VGPe[i][t] - VNaG);
			ITG = gTG * ainfG[i]*ainfG[i]*ainfG[i]*rG[i][t] * (VGPe[i][t] - VCaG);
			ICaG = gCaG * sinfG[i]*sinfG[i] * (VGPe[i][t] - VCaG);
			IAHPG = gAHPG * (VGPe[i][t] - VKG) * (CaG[i][t] / (CaG[i][t] + k1G));

			// GPe --> GPe synaptic current
			for (int k = STN_Num; k < STN_Num + GPe_Num; k++){
				if (adj_matrix[k][i + STN_Num] == -1){
					IsynGG = IsynGG + (strength_matrix[k][i + STN_Num] * sG[k - STN_Num][t] * (VGPe[i][t] - VsynGG));
					counterGG++;
				}
			}

			// STN --> GPe synaptic current
			for (int j = 0; j < STN_Num; j++){
				if (adj_matrix[j][i + STN_Num] == 1){
    					IsynSG = IsynSG + (strength_matrix[j][i + STN_Num] * sS[j][t] * (VGPe[i][t] - VsynSG));
					counterSG++;
				}
			}

			IsynG = (IsynSG / float(counterSG)) + (IsynGG / float(counterGG));
			LFP_GPe = LFP_GPe + VGPe[i][t];

			IsynSG = 0.0;
			IsynGG = 0.0;
			counterGG = 0;
			counterSG = 0;

			VGPe[i][t + 1] = VGPe[i][t] + (dt / Cm) * (- ILG - IKG - INaG - ITG - ICaG - IAHPG - IsynG + IappG[i] + I_stim_GPe);
			sG[i][t + 1] = sG[i][t] + dt * (alphaG * HinfG[i] * (1.0 - sG[i][t]) - betaG * sG[i][t]);
		        nG[i][t + 1] = nG[i][t] + dt * (phinG * (ninfG[i] - nG[i][t]) / taunG[i]);
			hG[i][t + 1] = hG[i][t] + dt * (phihG * (hinfG[i] - hG[i][t]) / tauhG[i]);
			rG[i][t + 1] = rG[i][t] + dt * (phirG * (rinfG[i] - rG[i][t]) / taurG);
			CaG[i][t + 1] = CaG[i][t] + dt * (epsG * (-ICaG - ITG - kCaG * CaG[i][t]));

			if (VGPe[i][t] > 0.0){
				if (VGPe[i][t] > VGPe[i][t - 1] && VGPe[i][t] > VGPe[i][t + 1]){
					tf_GPe[i] = t;
					outfile3 << t * dt / 1000. << '\t' << i + 1 << endl;
            			}  
			}
		}

//**************************************************************************************************

		// STN Cell
		for (int i = 0; i < STN_Num; i++){
			ainfS[i] = 1.0 / (1.0 + exp(-(VSTN[i][t] - thaS) / sigaS));
			sinfS[i] = 1.0 / (1.0 + exp(-(VSTN[i][t] - thsS) / sigsS));
			rinfS[i] = 1.0 / (1.0 + exp(-(VSTN[i][t] - thrS) / sigrS));
			minfS[i] = 1.0 / (1.0 + exp(-(VSTN[i][t] - thmS) / sigmS));
			ninfS[i] = 1.0 / (1.0 + exp(-(VSTN[i][t] - thnS) / signS));
			hinfS[i] = 1.0 / (1.0 + exp(-(VSTN[i][t] - thhS) / sighS));
			HinfS[i] = 1.0 / (1.0 + exp(-(VSTN[i][t] - thS - thHS) / sigHS));
			binfS[i] = 1.0 / (1.0 + exp((rS[i][t] - thetab) / sigmab)) - 1.0 / (1.0 + exp(-thetab / sigmab));	
			taunS[i] = taun0S + taun1S / (1.0 + exp(-(VSTN[i][t] - thnSt) / snS));
			tauhS[i] = tauh0S + tauh1S / (1.0 + exp(-(VSTN[i][t] - thhSt) / shS));
			taurS[i] = taur0S + taur1S / (1.0 + exp(-(VSTN[i][t] - thrSt) / srS));

			ILS = gLS * (VSTN[i][t] - VLS);
			IKS = gKS * nS[i][t]*nS[i][t]*nS[i][t]*nS[i][t] * (VSTN[i][t] - VKS);
			INaS = gNaS * minfS[i]*minfS[i]*minfS[i]*hS[i][t] * (VSTN[i][t] - VNaS);
			ITS = gTS * ainfS[i]*ainfS[i]*ainfS[i]*binfS[i]*binfS[i] * (VSTN[i][t] - VCaS);
			ICaS = gCaS * sinfS[i]*sinfS[i] * (VSTN[i][t] - VCaS);
			IAHPS = gAHPS * (VSTN[i][t] - VKS) * (CaS[i][t] / (CaS[i][t] + k1S));

			// GPe --> STN synaptic current
			for (int j = STN_Num; j < STN_Num + GPe_Num; j++){
				if (adj_matrix[j][i + STN_Num] == -1){
    					IsynGS = IsynGS + (strength_matrix[j][i + STN_Num] * sG[j - STN_Num][t] * (VSTN[i][t] - VsynGS));
					counterGS++;
				}
			}

			IsynS = (IsynGS / float(counterGS));
			LFP_STN = LFP_STN + VSTN[i][t];

			IsynGS = 0.0;
			counterGS = 0;

			VSTN[i][t + 1] = VSTN[i][t] + (dt / Cm) * (- ILS - IKS - INaS - ITS - ICaS - IAHPS - IsynS + IappS[i] + I_stim_STN);
			sS[i][t + 1] = sS[i][t] + dt * (alphaS * HinfS[i] * (1.0 - sS[i][t]) - betaS * sS[i][t]);
			nS[i][t + 1] = nS[i][t] + dt * (phinS * (ninfS[i] - nS[i][t]) / taunS[i]);
			hS[i][t + 1] = hS[i][t] + dt * (phihS * (hinfS[i] - hS[i][t]) / tauhS[i]);
			rS[i][t + 1] = rS[i][t] + dt * (phirS * (rinfS[i] - rS[i][t]) / taurS[i]);
			CaS[i][t + 1] = CaS[i][t] + dt * (epsS * (-ICaS - ITS - kCaS * CaS[i][t]));

			if (VSTN[i][t] > 0.0){
				if (VSTN[i][t] > VSTN[i][t - 1] && VSTN[i][t] > VSTN[i][t + 1]){
					tf_STN[i] = t;
					outfile4 << t * dt / 1000. << '\t' << i + 1 << endl;
            			} 
        		}
		}

//*************************************************************************************************************

	// iSTDP
	for (int i = STN_Num; i < STN_Num + GPe_Num; i++){
		for (int j = 0; j < STN_Num; j++){
			if (t == tf_GPe[i - STN_Num] && adj_matrix[i][j] == -1){
				delta_t = abs(tf_GPe[i - STN_Num] - tf_STN[j]) * dt;
				delta_g = eta_STDP * (exp(- delta_t / tau_STDP) - alpha_STDP);
				strength_matrix[i][j] = strength_matrix[i][j] + delta_g;
			}

	if (strength_matrix[i][j] > gmax) strength_matrix[i][j] = gmax;
	if (strength_matrix[i][j] < gmin) strength_matrix[i][j] = gmin;

		}
	}

//*************************************************************************************************************
		
		outfile5 << t * dt / 1000. << '\t' <<  LFP_GPe / float (GPe_Num) << endl;
		outfile6 << t * dt / 1000. << '\t' << -LFP_STN / float (STN_Num) << endl;

		LFP_GPe = 0.0;
		LFP_STN = 0.0;

	}// End of the time loop.	

//*************************************************************************************************************
//*************************************************************************************************************

	outfile1.close();
	outfile2.close();
	outfile3.close();
	outfile4.close();
	outfile5.close();
	outfile6.close();

return 0;
}

//*************************************************************************************************************
//*************************************************************************************************************
