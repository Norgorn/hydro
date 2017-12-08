package ru.norgorn.fiz1;

import static java.lang.Math.pow;

import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.function.BiConsumer;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class Fz2Main implements Runnable{
	
	int stepsT = 40_000;
	int stepsX = 50;
	int stepsZ = 50;
	int lastTInd = stepsT-1;
	int lastXInd = stepsX-1;
	int lastZInd = stepsZ-1;
	double maxT = 0.5;
	double L = 1;
	double maxX = L;
	double maxZ = 0.3;
	double dt = maxT/stepsT;
	double dx = maxX/stepsX;
	double dz = maxZ/stepsZ;
	
	double C = 0.3;
	double Q = 0.8;
	double F = 0.4;
	double A0 = 1;
	double P = 0.1;
	double D = 1;
	double p = 2;
	double n = 1;
	double g = 1;
	double beta_p = 1;
	double alpha = 0.5;
	double beta = 0.5;
	
	
	Fz2Values currentValues;
	Fz2Values previousValues;
	
	DecimalFormat format = new DecimalFormat("#0.0000");
	PrintWriter cOut;
	PrintWriter qOut;
	PrintWriter kOut;

	public static void main(String[] args) {
		try{
			new Fz2Main().run();	
		} catch(Exception e){
			e.printStackTrace();
		}
	}
	
	@Override
	public void run(){
		try{
			String cPath = cSavePath();
			String qPath = qSavePath();
			String kPath = kSavePath();
			new PrintWriter(cPath).close();
			new PrintWriter(qPath).close();
			cOut = new PrintWriter(cPath);
			qOut = new PrintWriter(qPath);
			kOut = new PrintWriter(kPath);
		} catch (Exception e){
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

		iterateT((t,i)->{
			if( i> 1 && i % (stepsT/1000) == 0 ){
				printRow(previousValues.c[10], i);
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
			iterateX((x,j) -> {
				iterateZ((z,m)->{
					next_v(j, m, z);
				});
			});
			iterateX((x,j)->{
				iterateZ((z,m)->{
					if(j==0)
						currentValues.c[j][m] = leftBoundaryState(previousValues.c, z, m);
					else if ( j==lastXInd)
						currentValues.c[j][m] = rightBoundaryState(previousValues.c, z, m);
					else{
						double r = next_c(j, m);
						currentValues.c[j][m] = r;
					}
				});
			});
			
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
			return C;
		return 0;
	}

	double leftBoundaryState(double[][] c, double z, int m){
		return C;
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
		double kk = A0*pow(F-qq,3) / pow(F-qq-1,2);
		return kk;
	}
	
	protected void next_v(int j, int m, double z) {
		
		double[][] k = previousValues.k;
		double[][] c = previousValues.c;
		
		double firctDerCX  = firstDerX(c, j, m);
		double vx = -k[j][m]/n*(P - p*g*z*beta_p*firctDerCX);
		//double vx = -k[j][m]/n*P;
		currentValues.vx[j][m] = vx;
		
		if(m ==0 || m==lastZInd){
			currentValues.vz[j][m] = 0;
		}
		else{
			double firctDerCZ  = firstDerZ(c, j, m);
			double vz = -k[j][m]/n* p*g* ( 1 - beta_p*c[j][m] - beta_p*z*firctDerCZ);

//			double vz = -k[j][m]/n * p*g * (1 - beta_p*c[j][m]);
			
			currentValues.vz[j][m] = vz;
		}
	}
	
	protected double next_q(int j, int m) {
		double[][] q = previousValues.q;
		double[][] c = previousValues.c;
		
		double qq =  q[j][m]+dt*(alpha*(Q-q[j][m])*c[j][m]-beta*q[j][m]);
		return qq;
	}
	
	protected double next_c(int j, int m) {
		double[][] q = previousValues.q;
		double[][] nextQ = currentValues.q;
		double[][] c = previousValues.c;
		double[][] vx = previousValues.vx;
		double[][] vz = previousValues.vz;
		
		double firctDerCX  = firstDerX(c, j, m);
		double firctDerCZ  = firstDerZ(c, j, m);
		double k1 = q[j][m]-nextQ[j][m];
		double k2 = dt*D*(secondDerX(c,j,m)+secondDerZ(c,j,m));
		double k3 = - dt*(vx[j][m]*firctDerCX + vz[j][m]*firctDerCZ);
		double kc = c[j][m];
		double r = kc + k1+ k2 + k3;
		return r;
	}
	
	protected double firstDerX(double[][] c, int j, int m) {
		if( j ==0)
			return 0;
		if (j == lastXInd)
			return c[j-1][m];
		return (c[j+1][m]-c[j-1][m])/(2*dx); 
	}
	
	protected double secondDerX(double[][] c, int j, int m) {
		double r = (c[j+1][m]-2*c[j][m]+c[j-1][m]) / pow(dx, 2);
		return r;
	}
	
	protected double firstDerZ(double[][] c, int j, int m) {
		if(m==0)
			return 0;
		if (m == lastZInd)
			return 0;
		return (c[j][m+1]-c[j][m-1])/(2*dz); 
	}
	
	protected double secondDerZ(double[][] c, int j, int m) {
		if(m==0)
			return 0;
		if (m == lastZInd)
			return 0;
		double r = (c[j][m+1]-2*c[j][m]+c[j][m-1]) / pow(dz, 2);
		return r;
	}
	
	protected String cSavePath() {
		return "C:\\Users\\Sunny\\Documents\\Wolfram Mathematica\\2d\\c.txt";
	}
	
	protected String qSavePath() {
		return "C:\\Users\\Sunny\\Documents\\Wolfram Mathematica\\2d\\q.txt";
	}
	
	protected String kSavePath() {
		return "C:\\Users\\Sunny\\Documents\\Wolfram Mathematica\\2d\\k.txt";
	}
	
	void iterateT(BiConsumer<Double, Integer> operation){
		double t=0;
		for(int i=0; i<=lastTInd; i++){
			operation.accept(t, i);
			t+=dt;
		}
	}
	
	void iterateX(BiConsumer<Double, Integer> operation){
		double x=0;
		for(int i=0; i<= lastXInd; i++){
			operation.accept(x, i);
			x+=dx;
		}
	}
	
	void iterateZ(BiConsumer<Double, Integer> operation){
		IntStream.range(0, stepsZ).parallel()
			.forEach(m -> {
				double z = m*dz;
				operation.accept(z, m);
			});
	}
	
	protected void printRow(double[] c, int i) {
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
	
	public String formatRow(double[] c) {
		String row = Arrays.stream(c).boxed()
				.map(d ->format.format(d)).collect(Collectors.joining("\t"))
				.replaceAll(",", "\\.");
		return row;
	}
}
