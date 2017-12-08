package ru.norgorn.fiz1;

import static java.lang.Math.pow;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class Fz1MainNonDimensional extends Fz1Main {
	
	double Q = 1.5;
	double Pe = 1;
	double F = 1;
	double a = 0.5;
	double b = 0.5;
	
	double[] k;
	
	static String params;
	static List<Double> stoppedAt = new ArrayList<>();
	
	public static void main(String[] args) {
		//double[] QQ = new double[] {1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 1.95, 2, 2.05, 2.1, 2.15, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8};
		double[] QQ = new double[] {8,9,10, 11, 12, 13 ,14 ,15 ,16 ,17};
		//double[] aa = new double[] {0.9, 1, 10};
		for(int i=0; i< QQ.length; i++){
			Fz1MainNonDimensional p = new Fz1MainNonDimensional();
			p.a = QQ[i];
			p.run();
		}
		System.out.println(params);
		AtomicInteger ind = new AtomicInteger(0);
		String str = stoppedAt.stream()
			.map(time -> "{"+QQ[ind.getAndIncrement()]+","+time+"}")
			.collect(Collectors.joining(","));
		System.out.println("Stopped at: "+str);
	}
	
	@Override
	public void run(){
		params = "a="+a+"_b="+b+"_Q="+Q+"_F="+F;
		
		k = new double[stepsT];
		super.run();
		
		try{
			String kPath = "C:\\Users\\Sunny\\Documents\\Wolfram Mathematica\\k_nd_"+params+".txt";
			new PrintWriter(kPath).close();
			PrintWriter kOut = new PrintWriter(kPath);
			double[] kReduced = IntStream.range(0, k.length).filter(i -> i % (stepsT/500) == 0)
				.mapToDouble(i -> k[i]).toArray();
			printRow(kReduced);
			saveRow(kReduced, kOut);
			kOut.close();
		} catch(Exception e){
			e.printStackTrace();
		}
	}
	
	@Override
	protected void timeStepFinished(int i, double t) {
		double qq = q[i][lastXInd];
		k[i] = pow(F-qq,3) / pow(F-qq-1/C,2);
		if(!doContinue(t, i)){
			System.out.println("STOPPED AT T: "+t+" (step "+i+")");
			stoppedAt.add(t);
		} else if(i == lastTInd-1) {
			stoppedAt.add(100.0);
		}
	}
	
	@Override
	double initialState(double x){
		return 0;
	}

	@Override
	double boundaryState(boolean leftBoundary, double[] c){
		if(leftBoundary)
			return 1;
		else
			return c[lastXInd-1];
	}
	
	@Override
	protected double next_q(double[][] c, double[][] q, Integer i, Integer j) {
		double qq =  q[i][j]+dt*(a*(Q-q[i][j])*c[i][j]-b*q[i][j]);
		return qq;
	}
	
	@Override
	protected double next_c(double[][] c, double[][] q, Integer i, Integer j) {
		double k1 = q[i][j]-q[i+1][j];
		double k2 = dt*secondDer(c,i,j);
		double k3 = Pe*pow(F-q[i][j],3)/pow(F-1/C-q[i][j],2)*dt*firstDer(c,i,j);
		double kc = c[i][j];
		double r = kc + k1+ k2 + k3;
		return r;
	}
	
	@Override
	protected boolean doContinue(double t, int i){
		return k[i]/k[0] > 0.1;
		//return true;
	}
	
	@Override
	protected String cSavePath() {
		return "C:\\Users\\Sunny\\Documents\\Wolfram Mathematica\\c_nd_"+params+".txt";
	}
	
	@Override
	protected String qSavePath() {
		return "C:\\Users\\Sunny\\Documents\\Wolfram Mathematica\\q_nd_"+params+".txt";
	}
}
