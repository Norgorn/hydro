package ru.norgorn.fiz1;

public class Fz1Main2 {
	
	int stepsT = 100;
	int stepsX = 10;
	int lastTInd = stepsT-1;
	int lastXInd = stepsX-1;
	double maxT = 1;
	double L = 1;
	double maxX = L;
	double dt = maxT/stepsT;
	double dx = maxX/stepsX;
	double D =0.2;
	double C =0;
	
	public static void main(String[] args){
		try{
			new Fz1Main2().run();
		} catch (Exception e){
			e.printStackTrace();
		}
	}
	
	public void run(){
		double[][] c = new double[stepsT][stepsX];
		for(int j=1; j< stepsX; j++){
			c[0][j] = 0;
		}
		c[0][5] = 10;
		for(int i=0; i< stepsT-1; i++){
			System.out.print(i+1+" ");
			for(int j=0; j< stepsX; j++){
				double r;
				if(j==0)
					r = C;
				else if (j==lastXInd)
					r = c[i][j-1];
				else{
					double r1 = (c[i][j+1]-2*c[i][j]+c[i][j-1]);
					r = c[i][j] + D*r1;//*dt/Math.pow(dx, 2);
				}
				c[i+1][j] = r;
				System.out.print(r+" ");
			}
			System.out.println();
		}
	}
}
