package ru.norgorn.fiz1;

import static java.lang.Math.pow;

import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.stream.Collectors;

public class Fz2MainNonDimensional extends FzDerivatives implements Runnable{
	
	public static double C0 = 0.2;
	public static double Pe = 100;
	public static double Fi = 1;
	public static double Q = 1.5;
	public static double Rp = 100;
	public static double beta = 0.1d;
	public static double a = 50;
	public static double b = 0;
	public static double L = 1;
	
	
	Fz2Values currentValues;
	Fz2Values previousValues;
	FzPoisson2D poisson = new FzPoisson2D();
	
	static DecimalFormat formatGeneral = new DecimalFormat("#0.000");
	static DecimalFormat formatPrecise = new DecimalFormat("#0.00000000");
	PrintWriter cOut;
	PrintWriter qOut;
	PrintWriter kOut;
	PrintWriter vxOut;
	PrintWriter vzOut;
	PrintWriter backflowOut;
	static long start = System.nanoTime();

	public static void main(String[] args) {
		try{
			new Fz2MainNonDimensional().run();
		} catch(Exception e){
			e.printStackTrace();
		}
	}
	
	@Override
	public void run(){
		maxT = 0.2;
		stepsT = 100_000;
		stepsX = 50;
		stepsZ = 50;
		init();
		
		try{
			String cPath = cSavePath();
			String qPath = qSavePath();
			String kPath = kSavePath();
			String vxPath = vxSavePath();
			String vzPath = vzSavePath();
			String backPath = backSavePath();
			new PrintWriter(cPath).close();
			new PrintWriter(qPath).close();
			System.out.println(cPath);
			cOut = new PrintWriter(cPath);
			qOut = new PrintWriter(qPath);
			kOut = new PrintWriter(kPath);
			vxOut = new PrintWriter(vxPath);
			vzOut = new PrintWriter(vzPath);
			backflowOut = new PrintWriter(backPath);
		} catch (Throwable e){
			e.printStackTrace();
		}
		
		
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
				printRow(previousValues.c[point], i / stepNumToShow);
//				printRow(previousValues.backflow[point], i / stepNumToShow, formatPrecise);
//				printRow(previousValues.q[point], i / stepNumToShow);
//				printRow(previousValues.p[point], i / stepNumToShow);
//				printRow(previousValues.vx[point], i / stepNumToShow);
//				printRow(previousValues.vz[point], i / stepNumToShow);
//				printRow(previousValues.k[point], i / stepNumToShow);
				saveRow(previousValues.c, cOut);
				saveRow(previousValues.q, qOut);
				saveRow(previousValues.k, kOut);
				saveRow(previousValues.vx, vxOut);
				saveRow(previousValues.vz, vzOut);
				saveRow(previousValues.backflow, backflowOut, formatPrecise);
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
		FzPoisson2DExplicit poisson = new FzPoisson2DExplicit(previousValues, currentValues, Pe, Rp, C0, beta);
		poisson.maxT = maxT;
		poisson.stepsT = stepsT;
		poisson.stepsX = stepsX;
		poisson.stepsZ = stepsZ;
		poisson.run();
		return poisson.currentValues.p;
	}
	
	protected void next_v(int j, int m, double x, double z) {
		double[][] kk = previousValues.k;
		double[][] cc = previousValues.c;
		double[][] p = previousValues.p;
		
		double k = kk[j][m];
		double c = cc[j][m];
		
		double vx;
		if(j == 0)
			vx = + Pe;
		else
			vx = - Pe*k*(firstDerX(p, j, m)) + Rp*k*c;
		//double vx = Pe*k*(1 - firstDerX(p, j, m)) + Rp*k*z*firstDerX(cc, j, m); // Complicated version, excluded now
		currentValues.vx[j][m] = vx;
		if(vx <0)
			getClass();
		
		if(m ==0 || m==lastZInd){
			currentValues.vz[j][m] = 0;
		}
		else{
			double vz = - Pe*k*firstDerZ(p, j, m) + Rp*k*c; 
			//double vz =  - k*Pe * firstDerZ(p, j, m) + Rp*k*c + Rp*k*z*firstDerZ(cc, j, m); // Complicated version, excluded now
			//vx = vx - Rp*k/beta/C0;   // EXCLUDED FROM EQUATIONS
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
	
	protected String cSavePath() {
		return "C:\\Users\\Sunny\\Documents\\Wolfram Mathematica\\2d_nd\\c_"+parametersString()+".txt";
	}

	protected String qSavePath() {
		return "C:\\Users\\Sunny\\Documents\\Wolfram Mathematica\\2d_nd\\q_"+parametersString()+".txt";
	}
	
	protected String kSavePath() {
		return "C:\\Users\\Sunny\\Documents\\Wolfram Mathematica\\2d_nd\\k_"+parametersString()+".txt";
	}
	
	protected String backSavePath() {
		return "C:\\Users\\Sunny\\Documents\\Wolfram Mathematica\\2d_nd\\backflow_"+parametersString()+".txt";
	}
	
	protected String vSavePath() {
		return "C:\\Users\\Sunny\\Documents\\Wolfram Mathematica\\2d_nd\\v_"+parametersString()+".txt";
	}
	
	protected String vxSavePath() {
		return "C:\\Users\\Sunny\\Documents\\Wolfram Mathematica\\2d_nd\\vx_"+parametersString()+".txt";
	}
	
	protected String vzSavePath() {
		return "C:\\Users\\Sunny\\Documents\\Wolfram Mathematica\\2d_nd\\vz_"+parametersString()+".txt";
	}
	
	private String parametersString() {
		return "C0="+C0+"_Fi="+Fi+"_Q="+Q+"_a="+a+"_b="+b+"_Pe="+Pe+"_Rp="+Rp+"_beta="+beta;
	}
	
	protected static void printRow(double[] c, int i) {
		printRow(c, i, formatGeneral);
	}
	
	protected static void printRow(double[] c, int i, DecimalFormat format) {
		String row = formatRow(c, format);
		long end = System.nanoTime();
		long time = (end - start)/1_000_000;
		start = end;
		System.out.println(i+"("+time+"): "+row);
	}
	
	protected void saveRow(double[][] row, PrintWriter out) {
		saveRow(row, out, formatGeneral);
	}
	
	protected void saveVectorRow(double[][] vx, double[][] vz, PrintWriter out) {
		String[][] matrix = new String[stepsX][stepsZ];
		iterateX((x,j)-> {
			iterateZ((z,m)-> {
				matrix[j][m] = "{"+formatGeneral.format(vx[j][m]).replaceAll(",", "\\.")
						+","+formatGeneral.format(vz[j][m]).replaceAll(",", "\\.")+"}"; 
			});
		});
		
		String rowStr = Arrays.stream(matrix)
			.map(r ->Arrays.stream(r).collect(Collectors.joining(",")))
			.collect(Collectors.joining(","));
		
		try {
			out.println(rowStr);
			out.flush();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	
	protected void saveRow(double[][] row, PrintWriter out, DecimalFormat format) {
		StringBuilder builder = new StringBuilder();
		for (int i=0; i< row.length; i++){
			builder.append(formatRow(row[i], format)).append(" ");
		}
		
		try {
			out.println(builder.toString());
			out.flush();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static String formatRow(double[] c, DecimalFormat format) {
		String row = Arrays.stream(c).boxed()
				.map(d ->format.format(d)).collect(Collectors.joining("\t"))
				.replaceAll(",", "\\.");
		return row;
	}
	
//	private double[][] next_p() {
//	poisson = new FzPoisson2D();
//	poisson.stepsX = stepsX;
//	poisson.stepsZ = stepsZ;
//	double[][] k = currentValues.k;
//	double[][] c = currentValues.c;
//	double[][] p = currentValues.p;
//	double[] pvc = previousValues.c[2];
//	
//	poisson.firstDerCoefFunctionX = (x,j,z,m) -> {
//		double val = 1/k[j][m] * firstDerX(k, j, m);
//		return val;
//	};
//	poisson.firstDerCoefFunctionZ = (x,j,z,m) -> {
//		double val = 1/k[j][m] * firstDerZ(k, j, m);
//		return val;
//	};
//	poisson.rightDependentKoefFunctionX = (x,j,z,m) -> {
//		double val = 1/k[j][m] * firstDerZ(k, j, m)*firstDerZ(p, j, m);
//		return val;
//	};
//	poisson.rightDependentKoefFunctionZ = (x,j,z,m) -> {
//		double val = 1/k[j][m] * firstDerX(k, j, m)*firstDerX(p, j, m);
//		return val;
//	};
//	
//	poisson.rightSideFunction = (x,j,z,m) -> {
//		double kk = k[j][m];
//		double cc = c[j][m];
//		double dkx = firstDerX(k, j, m);
//		double dkz = firstDerZ(k, j, m);
//		double dcx = firstDerX(c, j, m);
//		double dcz = firstDerZ(c, j, m);
////		double val = Rp/Pe * (kk*z*secondDerX(c, j, m) + 2*kk*dcz + kk*z*secondDerZ(c, j, m))
////				+  -dkx + Rp*dkx*z*firstDerX(c, j, m)
////				+ -Rp/Pe/(beta*C0)*dkz + Rp/Pe*dkz*c[j][m] + Rp/Pe*dkz*z*dcz;
//		double val = - 1/kk*dkx
//				+ Rp/Pe * (dkx/kk*z*dcx + z*secondDerX(c, j, m) + dkz/kk*(cc + z*dcz - 1/beta/C0) + 2*dcz + z*secondDerZ(c, j, m));
//		return val;
//	};
//	
//	poisson.run();
//	return poisson.Pc;
//}
}
