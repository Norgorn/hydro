package ru.norgorn.fiz1;

import static java.lang.Math.pow;

import java.text.DecimalFormat;

public class Fz2MainNonDimensional extends FzDerivatives implements Runnable{
	
	private static final String basePath = "C:\\Users\\Sunny\\Documents\\Wolfram Mathematica\\2d_nd\\";
	//private static final String basePath = "C:\\Users\\Sunny\\Documents\\Wolfram Mathematica\\2d_nd_iteration_2\\";
	
	public static DecimalFormat formatPrecise = new DecimalFormat("#0.00000000");
	
	public static double C0 = 0.2;
	public static double Fi = 3;
	public static double Q = 5;
	public static double Rp = 25;
	public static double Pe = 1000;
	public static double a = 15;
	public static double b = 19;
	public static double L = 1;
	
	
	Fz2Values currentValues;
	Fz2Values previousValues;
	
	Fz2Printer cOut;
	Fz2Printer tOut;
	Fz2Printer qOut;
	Fz2Printer kOut;
	Fz2Printer vxOut;
	Fz2Printer vzOut;
	Fz2Printer backflowOut;

	public static void main(String[] args) {
		try{
			double[] rps = new double[]{10, 15, 20, 25, 30, 40, 45, 50, 55};
			//double[] rps = new double[]{30};
			//double[] pes = new double[]{10, 50, 100, 150, 200, 275, 350, 500, 750, 1000};
			double[] pes = new double[]{100};
			for(double rp : rps){
				Rp = rp;
				for(double pe : pes){
					Pe = pe;
					new Fz2MainNonDimensional().run();	
				}
			}
		} catch(Exception e){
			e.printStackTrace();
		}
	}
	
	@Override
	public void run(){
		maxT = 0.1;
		stepsT = 50_000;
		stepsX = 50;
		stepsZ = 50;
		init();
		
		cOut = Fz2Printer.clearFileAndPrepareWriter(basePath, parametersString(), "c");
		qOut = Fz2Printer.clearFileAndPrepareWriter(basePath, parametersString(), "q");
		kOut = Fz2Printer.clearFileAndPrepareWriter(basePath, parametersString(), "k");
		vxOut = Fz2Printer.clearFileAndPrepareWriter(basePath, parametersString(), "vx");
		vzOut = Fz2Printer.clearFileAndPrepareWriter(basePath, parametersString(), "vz");
		backflowOut = Fz2Printer.clearFileAndPrepareWriter(basePath, parametersString(), "backflow");		
		
		previousValues = new Fz2Values(stepsX, stepsZ);
		iterateX((x,j) -> {
			iterateZ((z,m)-> {
				previousValues.c[j][m] = initialState(x,z);
				previousValues.q[j][m] = 0;
				previousValues.vx[j][m] = 0;
				previousValues.vz[j][m] = 0;
				previousValues.k[j][m] = next_k(j, m);
				previousValues.backflow[j][m] = 0;
			});
		});
		
		double kAtLeft = previousValues.k[0][0];

		int stepNumToShow = (stepsT/1000);
		//int stepNumToShow = (stepsT/stepsT);
		iterateT((t,i)->{
			if( i> 1 && i % stepNumToShow  == 0 ){
				int point = 10;
				cOut.printRow(previousValues.c[point], i / stepNumToShow);
				cOut.printRow(previousValues.t[point], i / stepNumToShow);
//				cOut.printRow(previousValues.backflow[point], i / stepNumToShow, formatPrecise);
//				cOut.printRow(previousValues.q[point], i / stepNumToShow);
//				cOut.printRow(previousValues.p[point], i / stepNumToShow);
//				cOut.printRow(previousValues.vx[point], i / stepNumToShow);
//				cOut.printRow(previousValues.vz[point], i / stepNumToShow);
//				cOut.printRow(previousValues.k[point], i / stepNumToShow);
				cOut.saveRow(previousValues.c);
				qOut.saveRow(previousValues.q);
				kOut.saveRow(previousValues.k);
				vxOut.saveRow(previousValues.vx);
				vzOut.saveRow(previousValues.vz);
				backflowOut.saveRow(previousValues.backflow, formatPrecise);
			}
			
			currentValues = new Fz2Values(stepsX, stepsZ);

			if (i==lastTInd) {
				return;
			}

			iterateXParallel((x,j) -> {
				iterateZ((z,m)->{
					currentValues.q[j][m] = next_q(j, m);
				});
			});
			iterateXParallel((x,j) -> {
					iterateZ((z,m)->{
						if(j == -1){
							currentValues.k[j][m] = kAtLeft; 
						}
						else{
							currentValues.k[j][m] = next_k(j, m);
						}
					});
			});

			currentValues.p = next_p();
			
			iterateXParallel((x,j) -> {
				iterateZ((z,m)->{
					next_v(j, m, x, z);
				});
			});
			iterateXParallel((x,j)->{
				iterateZ((z,m)->{
					if(j==0)
						currentValues.c[j][m] = leftBoundaryState(previousValues.c, z, m);
					else if ( j==lastXInd)
						currentValues.c[j][m] = rightBoundaryState(previousValues.c, z, m);
					else if ( m==0){
						currentValues.c[j][m] = bottomBoundaryState(previousValues.c, z, j);
					}
					else if ( m==lastZInd){
						currentValues.c[j][m] = topBoundaryState(previousValues.c, z, j);
					}
					else{
						double r = next_c(j, m);
						currentValues.c[j][m] = r;
					}
				});
			});
			
			iterateXParallel((x,j)->{
				iterateZ((z,m)->{
					if(currentValues.c[j][m] < previousValues.c[j][m]){
						double backflowVal = previousValues.c[j][m] - currentValues.c[j][m];
						currentValues.backflow[j][m] = backflowVal;
					}
					if(m !=0){
						if(currentValues.c[j][m] > currentValues.c[j][m-1])
							currentValues.getClass();
					}
				});
			});
			previousValues = currentValues;
			timeStepFinished(i,t);
		});
	}

	protected void timeStepFinished(int i, double t) {
	}
	
	
	double initialState(double x, double z){
		if (x==0)
			return 1;
		return 0;
	}

	double leftBoundaryState(double[][] c, double z, int m){
		return 1;
	}
	
	double rightBoundaryState(double[][] c, double z, int m){
		return c[lastXInd-1][m];
	}
	
	private double topBoundaryState(double[][] c, double z, int j){
		return c[j][lastZInd-1];
	}
	
	private double bottomBoundaryState(double[][] c, double z, int j){
		return c[j][1];
	}
	
	protected double next_k(int j, int m) {
		double qq = previousValues.q[j][m];
		double kk = pow(Fi-qq,3) / pow(Fi-qq-1/C0, 2);
		return kk;
	}
	
	private double[][] next_p() {
		FzPoisson2DExplicit poisson = new FzPoisson2DExplicitPsi(previousValues, currentValues, Pe, Rp, C0);
		poisson.maxT = maxT;
		poisson.stepsT = stepsT;
		poisson.stepsX = stepsX;
		poisson.stepsZ = stepsZ;
		poisson.run();
		return poisson.currentValues.p;
	}
	
	protected void next_v(int j, int m, double x, double z) {
		double[][] p = previousValues.p;
		
		double vx;
		if(j == 0)
			vx = + Pe;
		else
			vx = firstDerZ(p, j, m);
		currentValues.vx[j][m] = vx;
		if(vx <0 )
			getClass();
		
		if(m ==0 || m==lastZInd){
			currentValues.vz[j][m] = 0;
		}
		else{
			//double vz = - Pe*k*firstDerZ(p, j, m) + Rp*k*c;
			double vz = - firstDerX(p, j, m); 
			currentValues.vz[j][m] = vz;
		}
	}
	
	protected double next_q(int j, int m) {
		double[][] q = previousValues.q;
		double[][] c = previousValues.c;
		
		double qq =  q[j][m]+dt*(a*(Q-q[j][m])*c[j][m]-b*q[j][m]);
		return qq;
	}
	
	protected double next_c(int j, int m) {
		double[][] q = previousValues.q;
		double[][] nextQ = currentValues.q;
		double[][] c = previousValues.c;
		double[][] vx = previousValues.vx;
		double[][] vz = previousValues.vz;
		
		double dcx  = firstDerX(c, j, m);
		double dcz  = firstDerZ(c, j, m);
		double k1 = - nextQ[j][m] + q[j][m];
		double k2 = dt*(secondDerX(c,j,m)+secondDerZ(c,j,m));
		double k3 = - dt*(vx[j][m]*dcx + vz[j][m]*dcz);
		double kc = c[j][m];
		double r = kc + k1+ k2 + k3;
		if(r <0)
			return 0;
		return r;
	}
	
	private String parametersString() {
		return "C0="+C0+"_Fi="+Fi+"_Q="+Q+"_a="+a+"_b="+b+"_Pe="+Pe+"_Rp="+Rp+"_beta=0.1";
	}
}
