package ru.norgorn.fiz1;

import static java.lang.Math.pow;

import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.stream.Collectors;

public class Fz2AlternateDirections extends FzDerivatives implements Runnable{
	
	public static double C0 = 0.2;
	public static double Pe = 2;
	public static double Fi = 2;
	public static double Q = 2;
	public static double Rp = 5;
	public static double beta = 0.1;
	public static double a = 0.5;
	public static double b = 0.5;
	
	
	Fz2Values currentValues;
	Fz2Values previousValues;
	FzPoisson2D poisson = new FzPoisson2D();
	
	static DecimalFormat format = new DecimalFormat("#0.000");
	PrintWriter cOut;
	PrintWriter qOut;
	PrintWriter kOut;

	public static void main(String[] args) {
		try{
			//new Fz2AlternateDirections().run();	
		} catch(Exception e){
			e.printStackTrace();
		}
	}
	
	@Override
	public void run(){
		maxT = 0.5;
		stepsT = 20_000;
		stepsX = 50;
		stepsZ = 50;
		init();
		
		try{
			String cPath = cSavePath();
			String qPath = qSavePath();
			String kPath = kSavePath();
			new PrintWriter(cPath).close();
			new PrintWriter(qPath).close();
			cOut = new PrintWriter(cPath);
			qOut = new PrintWriter(qPath);
			kOut = new PrintWriter(kPath);
		} catch (Throwable e){
			e.printStackTrace();
		}
		
		
		previousValues = new Fz2Values(stepsX, stepsZ);
		iterateX((x,i) -> {
			iterateZ((z,m)-> {
				previousValues.c[i][m] = initialState(x,z);
				previousValues.q[i][m] = 0;
				previousValues.vx[i][m] = 0;
				previousValues.vz[i][m] = 0;
			});
		});

		int stepNumToShow = (stepsT/1000);
		iterateT((t,i)->{
			if( i> 1 && i % stepNumToShow  == 0 ){
				printRow(previousValues.c[10], i / stepNumToShow);
				//printRow(previousValues.vx[10], i / stepNumToShow);
				//printRow(previousValues.vz[10], i / stepNumToShow);
				saveRow(previousValues.c, cOut);
				saveRow(previousValues.q, qOut);
				saveRow(previousValues.k, kOut);
			}
			
			currentValues = new Fz2Values(stepsX, stepsZ);

			if (i==lastTInd) {
				return;
			}

			iterateX((x,j) -> {
				iterateZ((z,m)->{
					currentValues.q[j][m] = next_q(j, m);
				});
			});
			iterateX((x,j) -> {
				iterateZ((z,m)->{
					currentValues.k[j][m] = next_k(j, m);
				});
			});

			currentValues.p = next_p();
			
			iterateX((x,j) -> {
				iterateZ((z,m)->{
					next_v(j, m, z);
				});
			});
			iterateXParallel((x,j)->{
				iterateZ((z,m)->{
					if(j==0)
						currentValues.c[j][m] = leftBoundaryState(previousValues.c, z, m);
					else if ( j==lastXInd)
						currentValues.c[j][m] = rightBoundaryState(previousValues.c, z, m);
				});
			});
			nextAllCValues();
			
			iterateX((x,j)->{
				iterateZ((z,m)->{
					if(currentValues.c[j][m] < previousValues.c[j][m] ||
							currentValues.q[j][m] < previousValues.q[j][m]){
						currentValues.getClass();
					}
					if(m !=0){
						if(currentValues.c[j][m] > currentValues.c[j][m-1] ||
								currentValues.q[j][m] > currentValues.q[j][m-1])
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
	
	double topBoundaryState(double[][] c, double z, int m){
		return c[lastXInd-1][m];
	}
	
	double bottomBoundaryState(double[][] c, double z, int m){
		return c[lastXInd-1][m];
	}
	
	protected double next_k(int j, int m) {
		double qq = previousValues.q[j][m];
		double kk = pow(Fi-qq,3) / pow(Fi-qq-1/C0,2);
		return kk;
	}
	
	private double[][] next_p() {
		poisson = new FzPoisson2D();
		poisson.stepsX = stepsX;
		poisson.stepsZ = stepsZ;
		double[][] k = previousValues.k;
		double[][] c = previousValues.c;
		
		poisson.firstDerCoefFunctionX = (x,j,z,m) -> {
			double val = firstDerX(k, j, m);
			return val;
		};
		poisson.firstDerCoefFunctionZ = (x,j,z,m) -> {
			double val = firstDerZ(k, j, m);
			return val;
		};
		poisson.rightSideFunction = (x,j,z,m) -> {
			
			//double val = Rp/Pe * (k*z*secondDerX(c, j, m) + k*firstDerZ(c, j, m) + k*z*secondDerZ(c, j, m));
			
			double kk = k[j][m];
			double dkx = firstDerX(k, j, m);
			double dkz = firstDerZ(k, j, m);
			double dcz = firstDerZ(c, j, m);
			double val = Rp/Pe * (kk*z*secondDerX(c, j, m) + 2*kk*dcz + kk*z*secondDerZ(c, j, m))
					+  -dkx + Rp*dkx*z*firstDerX(c, j, m)
					+ -Rp/Pe/(beta*C0)*dkz + Rp/Pe*dkz*c[j][m] + Rp/Pe*dkz*z*dcz;
			
			return val;
		};
		
		poisson.run();
		return poisson.Pc;
	}
	
	protected void next_v(int j, int m, double z) {
		
		double[][] k = previousValues.k;
		double[][] c = previousValues.c;
		double[][] p = previousValues.p;
		
		double firstDerCX  = firstDerX(c, j, m);
		double vx = -k[j][m] * (Pe + Pe*firstDerX(p, j, m) - C0*Rp*z*firstDerCX);
		currentValues.vx[j][m] = vx;
		
		
		
		if(m ==0 || m==lastZInd){
			currentValues.vz[j][m] = 0;
		}
		else{
			double firstDerCZ  = firstDerZ(c, j, m);
			double vz = -k[j][m]*Rp * ( 1/beta - C0*c[j][m] - C0*z*firstDerCZ) - k[j][m]*Pe * firstDerZ(p, j, m);
			
			currentValues.vz[j][m] = vz;
		}
	}
	
	protected double next_q(int j, int m) {
		double[][] q = previousValues.q;
		double[][] c = previousValues.c;
		
		double qq =  q[j][m]+dt*(a*(Q-q[j][m])*c[j][m]-b*q[j][m]);
		return qq;
	}
	
	protected void nextAllCValues() {
		double[][] tmpC = new double[stepsX][stepsZ];
		double[][] tmpD = new double[stepsX][stepsZ];
		double[][] vx = currentValues.vx;
		double[][] vz = currentValues.vz;
		double[][] pq = previousValues.q;
		double[][] cq = currentValues.q;
		{
			double[][] pc = previousValues.c;
			double b;
			b = -1 - dt/(dx*dx);
			iterateZParallel((z,m)->{
				double cc = dt/(2*dx*dx) - vx[0][m]*dt/(4*dx);
				double dd = -pc[0][m] - dt/2 * (pc[1][m] - (cq[0][m] - pq[0][m])/dt);
				tmpC[0][m] = cc/b;
				tmpD[0][m] = dd/b;
				iterateXExcludeBorders((x,j)->{
					double a,c,d;
					a = dt/(2*dx*dx) + vx[j][m]*dt/(4*dx);
					c = dt/(2*dx*dx) - vx[j][m]*dt/(4*dx);
					tmpC[j][m] = c / (b - a*tmpC[j-1][m]);

					double val;
					if(m != lastZInd && m!=0)
						val = vz[j][m] * firstDerZ(pc, j, m);//(pc[j][m+1]-pc[j][m-1]) /2/dz
					else
						val =0;
					d = -pc[j][m] - dt/2 * (secondDerZ(pc, j, m) - val) + (cq[j][m] - pq[j][m])/2;
					tmpD[j][m] = (d - a*tmpD[j-1][m]) / (b - a*tmpC[j-1][m]);
				});
			});
			iterateZParallel((z,m)->{
				double[][] td = tmpD;
				currentValues.c[lastXInd][m] = td[lastXInd][m];
				iterateBackwardXExcludeBorders((x,j)->{
					double v1 = tmpD[j][m];
					double v2 = tmpC[j][m]*currentValues.c[j+1][m];
					double v = v1 - v2;
					currentValues.c[j][m] = v;
				});
			});
		}
		
		{
			double[][] pc = currentValues.c;
			double b;
			b = -1 - dt/(dz*dz);
			iterateXExcludeBordersParallel((x,j)->{
				double cc = dt/(2*dz*dz) - vz[j][0]*dt/(4*dz);
				double dd = -pc[j][0] - dt/2 * (pc[j][1] - (cq[j][0] - pq[j][0])/dt);
				tmpC[j][0] = cc/b;
				tmpD[j][0] = dd/b;
				iterateZ((z,m)->{
					if(m==0)
						return;
					double a,c,d;
					a = dt/(2*dz*dz) + vz[j][m]*dt/(4*dz);
					c = dt/(2*dz*dz) - vz[j][m]*dt/(4*dz);
					tmpC[j][m] = c / (b - a*tmpC[j][m-1]);

					double val = vx[j][m]*firstDerX(pc, j, m);
					d = -pc[j][m] - dt/2 * (secondDerX(pc, j, m) - val) + (cq[j][m] - pq[j][m])/2;
					tmpD[j][m] = (d - a*tmpD[j][m-1]) / (b - a*tmpC[j][m-1]);
				});
			});
			iterateXExcludeBordersParallel((x,j)->{
				double[][] td = tmpD;
				currentValues.c[j][lastZInd] = td[j][lastZInd];
				iterateBackwardZ((z,m)->{
					if(m==lastZInd)
						return;
					double v1 = tmpD[j][m];
					double v2 = tmpC[j][m]*currentValues.c[j][m+1];
					double v = v1 - v2;
					currentValues.c[j][m] = v;
				});
			});
		}
	}
	
	protected String cSavePath() {
		return "C:\\Users\\Sunny\\Documents\\Wolfram Mathematica\\2d_nd_alt\\c.txt";
	}
	
	protected String qSavePath() {
		return "C:\\Users\\Sunny\\Documents\\Wolfram Mathematica\\2d_nd_alt\\q.txt";
	}
	
	protected String kSavePath() {
		return "C:\\Users\\Sunny\\Documents\\Wolfram Mathematica\\2d_nd_alt\\k.txt";
	}
	
	protected static void printRow(double[] c, int i) {
		String row = formatRow(c);
		System.out.println(i+": "+row);
	}
	
	protected void saveRow(double[][] row, PrintWriter out) {
		StringBuilder builder = new StringBuilder();
		for (int i=0; i< row.length; i++){
			builder.append(formatRow(row[i])).append(" ");
		}
		
		try {
			out.println(builder.toString());
			out.flush();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static String formatRow(double[] c) {
		String row = Arrays.stream(c).boxed()
				.map(d ->format.format(d)).collect(Collectors.joining("\t"))
				.replaceAll(",", "\\.");
		return row;
	}
}