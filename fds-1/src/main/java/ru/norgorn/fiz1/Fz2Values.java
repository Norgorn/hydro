package ru.norgorn.fiz1;

public class Fz2Values {
	public final double[][] k;
	public double[][] p;
	public final double[][] vx;
	public final double[][] vz;
	public final double[][] q;
	public final double[][] c;
	public final double[][] backflow;
	
	
	public Fz2Values(int stepsX, int stepsZ){
		k = new double[stepsX][stepsZ];
		p = new double[stepsX][stepsZ];
		vx = new double[stepsX][stepsZ];
		vz = new double[stepsX][stepsZ];
		q = new double[stepsX][stepsZ];
		c = new double[stepsX][stepsZ];
		backflow = new double[stepsX][stepsZ];
	}
}
