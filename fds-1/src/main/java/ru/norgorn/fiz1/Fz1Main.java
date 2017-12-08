package ru.norgorn.fiz1;

import static java.lang.Math.pow;

import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.function.BiConsumer;
import java.util.stream.Collectors;

public class Fz1Main implements Runnable {

	int stepsT = 200_000;
	int stepsX = 100;
	int lastTInd = stepsT-1;
	int lastXInd = stepsX-1;
	double maxT = 8;
	double L = 1;
	double maxX = L;
	double dt = maxT/stepsT;
	double dx = maxX/stepsX;
	double C = 0.1;
	double al = 1;
	double beta = 0.2;
	double Q0 = 0.1;
	double A0 = 1;
	double Fi0 = 0.3;
	double D = 1;
	double P = 1;
	double n = 10;
	double gamma = -A0/n*P*pow(Fi0,3)/pow(1-Fi0,2);
	double k = (Math.atan2(gamma/pow(1+gamma,2), -1/pow(1+gamma,2))+2*Math.PI)/L;

	DecimalFormat format = new DecimalFormat("#0.0000");
	PrintWriter cOut;
	PrintWriter qOut;
	
	double[][] c;
	double[][] q;

	public static void main(String[] args) {
		new Fz1Main().run();
	}

	@Override
	public void run() {
		try{
			String cPath = cSavePath();
			String qPath = qSavePath();
			new PrintWriter(cPath).close();
			new PrintWriter(qPath).close();
			cOut = new PrintWriter(cPath);
			qOut = new PrintWriter(qPath);

			c = new double[stepsT][];
			q = new double[stepsT][];
			iterateT((t,i) -> {
				c[i] = new double[stepsX];
				q[i] = new double[stepsX];
			});

			iterateX((x,i) -> {
				c[0][i] = initialState(x);
				q[0][i] = 0;
			});

			iterateTWithCheck((t,i)->{
				if(i % 1000 == 0){
					printRow(c[i]);
					saveRow(c[i],cOut);
					saveRow(q[i],qOut);
				}

				if (i==lastTInd) {
					return;
				}

				iterateX((x,j) -> {
					q[i+1][j] = next_q(c, q, i, j);
				});
				iterateX((x,j)->{
					if(j==0)
						c[i+1][j] = boundaryState(true, c[i]);
					else if ( j==lastXInd)
						c[i+1][j] = boundaryState(false, c[i]);
					else{
						double r = next_c(c, q, i, j);
						c[i+1][j] = r;
					}
				});
				timeStepFinished(i,t);
			});
			System.out.println("DONE!");
			cOut.close();
			qOut.close();
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
	protected boolean doContinue(double t, int i){
		return true;
	}
	
	protected void timeStepFinished(int i, double t) {
	}

	protected double next_q(double[][] c, double[][] q, Integer i, Integer j) {
		double qq =  q[i][j]+dt*(al*(Q0-q[i][j])*c[i][j]-beta*q[i][j]);
		return qq;
	}

	protected double next_c(double[][] c, double[][] q, Integer i, Integer j) {
		double k1 = q[i][j]-q[i+1][j];
		double k2 = dt*D*secondDer(c,i,j);
		double k3 = dt*P/n*A0*pow(Fi0-q[i][j],3)/pow(Fi0-1-q[i][j],2)*firstDer(c,i,j);
		double kc = c[i][j];
		double r = kc + k1+ k2 + k3;
		return r;
	}

	double initialState(double x){
		return 0;// C+Math.exp(gamma*x)*Math.sin(k*x);
	}

	double boundaryState(boolean leftBoundary, double[] c){
		if(leftBoundary)
			return C;
		else
			return c[lastXInd-1];
	}

	protected double secondDer(double[][] c, int i, int j) {
		double r = (c[i][j+1]-2*c[i][j]+c[i][j-1]) / pow(dx, 2);
		return r;
	}

	protected double firstDer(double[][] c, int i, int j) {
		return (c[i][j+1] - c[i][j-1]) / (2*dx); 
	} 

	void iterateT(BiConsumer<Double, Integer> operation){
		double t=0;
		for(int i=0; i<=lastTInd; i++){
			operation.accept(t, i);
			t+=dt;
		}
	}
	
	void iterateTWithCheck(BiConsumer<Double, Integer> operation){
		double t=0;
		for(int i=0; i<=lastTInd; i++){
			operation.accept(t, i);
			if(!doContinue(t,i))
				break;
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

	protected void printRow(double[] c) {
		String row = formatRow(c);
		//System.out.println(i+"\t"+row);
		System.out.println(row);
	}

	protected void saveRow(double[] row, PrintWriter out) {
		String rowStr = formatRow(row);
		try {
			out.println(rowStr);
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
	
	protected String cSavePath() {
		return "C:\\Users\\Sunny\\Documents\\Wolfram Mathematica\\c.txt";
	}
	
	protected String qSavePath() {
		return "C:\\Users\\Sunny\\Documents\\Wolfram Mathematica\\q.txt";
	}
}
