package ru.norgorn.fiz1;

import static  ru.norgorn.fiz1.Fz2MainNonDimensional.printRow;

import java.awt.RadialGradientPaint;
import java.util.function.BiFunction;

import javax.swing.JMenu;

public class FzPoisson2D extends FzDerivatives implements Runnable {

	private double tau = 0.000001;
	public double[][] Pc; //current
	private double[][] Pp; //previous
	private double sumError = 0;
	private double stopCriteria = 0.0000001;
	
	public QuFunction rightSideFunction;
	public QuFunction firstDerCoefFunctionX;
	public QuFunction firstDerCoefFunctionZ;
	public QuFunction rightDependentKoefFunctionX;
	public QuFunction rightDependentKoefFunctionZ;
	

	private int stepNum=0;

	public static void main(String[] args) {
		try {
			FzPoisson2D run = new FzPoisson2D();
			run.stepsX = 100;
			run.stepsZ = 100;
			run.run();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	@Override
	public void run() {
		
		init();
		
		long start = System.nanoTime();
		Pc = new double[stepsX][stepsZ];
		iterateX((x,i) -> {
			iterateZ((z,m) -> {
				if( x == 0)
					Pc[i][m] = leftBoubdaryCondition();
				else if (x == maxX)
					Pc[i][m] = rightBoubdaryCondition();
				else if (z == 0)
					Pc[i][m] = bottomBoubdaryCondition();
				else if (z == maxZ)
					Pc[i][m] = topBoubdaryCondition();
				else
					Pc[i][m] = initialState(x, z);
			});
		});

		double dx = dx();
		double dz = dz();
		do {
			Pp = Pc;
			Pc = new double[stepsX][stepsZ];
			iterateX((x,i) -> {
				iterateZ((z,m) -> {
					if( x == 0)
						Pc[i][m] = leftBoubdaryCondition();
					else if (x == maxX)
						Pc[i][m] = rightBoubdaryCondition();
					else if (z == 0)
						Pc[i][m] = bottomBoubdaryCondition();
					else if (z == maxZ)
						Pc[i][m] = topBoubdaryCondition();
				});
			});

			iterateZExcludeBordersParallel((z,m) -> {
				double F1 = - Pp[0][m] - tau/2 * secondDerZ(Pp, 0, m) - tau/2 * rightSide(0, 0, z, m);
				double[] c = new double[stepsZ];
				double[] f = new double[stepsZ];
				double Cxx =  tau/(2*dx*dx);
				double Bxx = -1 - tau/dx/dx;
				
				
				c[0] = Cxx/Bxx;
				f[0] = F1/Bxx;
				iterateXExcludeBorders((x,j) -> {
					double k =0;
					double r = 0;
					if(firstDerCoefFunctionX != null){
						k = firstDerCoefFunctionX.apply(x, j, z, m);
					}
					if(rightDependentKoefFunctionX != null){
						r = rightDependentKoefFunctionX.apply(x, j, z, m);
					}
					double Ax = tau/(2*dx *dx) - tau/(4*dx)*k;
					double Bx = -1 - tau/dx/dx;
					double Cx =  tau/(2*dx*dx) + tau/(4*dx)*k;
					
					double F = - Pp[j][m] - tau/2 * secondDerZ(Pp, j, m) - tau/2*(r + rightSide(x, j, z, m));
					double cc = c[j-1];
					double ff = f[j-1];
					c[j] = Cx / (Bx - Ax*cc);
					f[j] = (F - Ax*ff) / (Bx - Ax*cc);
				});

				for(int j =stepsX-2; j >0; j--){
					Pc[j][m] = f[j] - c[j]*Pc[j+1][m];
				}
			});

			iterateXExcludeBordersParallel((x,j) -> {
				double F1 = - Pc[j][0] - tau/2 * secondDerX(Pc, j, 0) - tau/2 * rightSide(x, j, 0, 0);
				double[] c = new double[stepsZ];
				double[] f = new double[stepsZ];
				double Czz = tau/(2*dz*dz);
				double Bzz = -1 - tau/dz/dz;
				
				c[0] = Czz/Bzz;
				f[0] = F1/Bzz;
				iterateZExcludeBorders((z,m) -> {
					double k =0;
					double r =0;
					if(firstDerCoefFunctionZ != null){
						k = firstDerCoefFunctionZ.apply(x, j, z, m);
					}
					if(rightDependentKoefFunctionX != null){
						r = rightDependentKoefFunctionX.apply(x, j, z, m);
					}
					double Az = tau/(2*dz*dz) - tau/(4*dz)*k;
					double Bz = -1 - tau/dz/dz;
					double Cz = tau/(2*dz*dz) + tau/(4*dz)*k;
					
					double F = - Pc[j][m] - tau/2 * secondDerX(Pc, j, m) - tau/2*(r + rightSide(x, j, z, m));
					double cc = c[m-1];
					double ff = f[m-1];
					c[m] = Cz / (Bz - Az*cc);
					f[m] = (F - Az*ff) / (Bz - Az*cc);
				});

				for(int m =stepsZ-2; m >0; m--){
					Pc[j][m] = f[m] - c[m]*Pc[j][m+1];
				}
			});
			
			sumError = countSquareSum();
		} while(keepGoing());
		//System.out.println( (System.nanoTime() - start)/1_000_000l);
	}

	private boolean keepGoing() {
		stepNum++;
		//printRow(Pc[20], stepNum);
		boolean ret = sumError  > stopCriteria;
		if(!ret){
			getClass();
		}
		return ret && stepNum < 300;
		//return true;
	}

	private double rightSide(double x, int j, double z, int m) {
		double val;
		if(rightSideFunction != null)
			val = rightSideFunction.apply(x, j, z, m);
		else
			val = - Math.exp(-(x*x + z*z));
		return val;
	}

	public double initialState(double x, double z){
		return 0;
	}

	public double leftBoubdaryCondition(){
		return 0;
	}

	public double rightBoubdaryCondition(){
		return 0;
	}

	public double topBoubdaryCondition(){
		return 0;
	}

	public double bottomBoubdaryCondition(){
		return 0;
	}
	
	private double countSquareSum() {
		double[] sum = new double[]{0};
		iterateX((x,j) -> {
			iterateZ( (z,m) ->{
				double value = Pp[j][m] - Pc[j][m];
				sum[0] += value*value;
			});
		});
		return sum[0];
	}
}


/*
  iterateZExcludeBorders((z,m) -> {
				double F1 = - Pc[0][m] - tau/2 * secondDerZ(Pc, 0, m) - tau/2 * rightSide(0, 0, z, m);
				double[] a = new double[Fz2MainNonDimensional.stepsX];
				double[] b = new double[Fz2MainNonDimensional.stepsX];
				a[1] = -Bx/Cx;
				b[1] = -F1/Cx;
				iterateXExcludeBorders((x,j) -> {
					double F = - Pc[j][m] - tau/2 * secondDerZ(Pc, j, m) - tau/2 * rightSide(x, j, z, m);
					F = -F;
					double aa = a[j];
					double bb = b[j];
					a[j+1] = -Bx/(Ax*aa+Cx);
					b[j+1] = (F - Ax*bb)/(Ax*aa+Cx);
				});

				for(int j =Fz2MainNonDimensional.stepsX-2; j >0; j--){
					Pc[j][m] = a[j+1] * Pc[j+1][m] + b[j+1];
				}
			});

			iterateXExcludeBorders((x,j) -> {
				double F1 = - Pc[j][0] - tau/2 * secondDerX(Pc, j, 0) - tau/2 * rightSide(x, j, 0, 0);
				double[] a = new double[Fz2MainNonDimensional.stepsZ];
				double[] b = new double[Fz2MainNonDimensional.stepsZ];
				a[1] = -Bz/Cz;
				b[1] = -F1/Cz;
				iterateZExcludeBorders((z,m) -> {
					double F = - Pc[j][m] - tau/2 * secondDerX(Pc, j, m) - tau/2 * rightSide(x, j, z, m);
					F = -F;
					double aa = a[m];
					double bb = b[m];
					a[m+1] = -Bz/(Az*aa+Cz);
					b[m+1] = (F - Az*bb)/(Az*aa+Cz);
				});

				for(int m =Fz2MainNonDimensional.stepsZ-2; m >0; m--){
					Pc[j][m] = a[m+1] * Pc[j][m+1] + b[m+1];
				}
			});




















			iterateZExcludeBorders((z,m) -> {
				double F1 = - Pc[0][m] - tau/2 * secondDerZ(Pc, 0, m) - tau/2 * rightSide(0, 0, z, m);
				double[] c = new double[Fz2MainNonDimensional.stepsZ];
				double[] f = new double[Fz2MainNonDimensional.stepsZ];
				c[0] = Cx;
				f[0] = F1;
				iterateXExcludeBorders((x,j) -> {
					double F = - Pc[j][m] - tau/2 * secondDerZ(Pc, j, m) - tau/2 * rightSide(x, j, z, m);
					double cc = c[j-1];
					double ff = f[j-1];
					c[j] = Cx - Ax*Bx/cc;
					f[j] = F - Ax*ff/cc;
				});

				for(int j =Fz2MainNonDimensional.stepsX-2; j >0; j--){
					Pc[j][m] = (f[j]-Bz*Pc[j+1][j])/c[j];
				}
			});

			iterateXExcludeBorders((x,j) -> {
				double F1 = - Pc[j][0] - tau/2 * secondDerX(Pc, j, 0) - tau/2 * rightSide(x, j, 0, 0);
				double[] c = new double[Fz2MainNonDimensional.stepsZ];
				double[] f = new double[Fz2MainNonDimensional.stepsZ];
				c[0] = Cz;
				f[0] = F1;
				iterateZExcludeBorders((z,m) -> {
					double F = - Pc[j][m] - tau/2 * secondDerX(Pc, j, m) - tau/2 * rightSide(x, j, z, m);
					double cc = c[m-1];
					double ff = f[m-1];
					c[m] = Cz - Az*Bz/cc;
					f[m] = F - Az*ff/cc;
				});

				for(int m =Fz2MainNonDimensional.stepsZ-2; m >0; m--){
					Pc[j][m] = (f[m]-Bz*Pc[j][m+1])/c[m];
				}
			});
 */