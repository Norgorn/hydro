package ru.norgorn.fiz1;

import java.util.Arrays;

import com.google.common.util.concurrent.AtomicDouble;

public class FzPoisson2DExplicit extends FzDerivatives implements Runnable  {
	
	private double stopCriteria = 0.0001;
	
	public final Fz2Values currentValues;
	private final Fz2Values previousValues;
	private final double Pe;
	private final double Rp;
	private final double C0;
	private final double beta;
	private final double L=1;
	
	private final boolean simpleCase = true;
	
	public FzPoisson2DExplicit(Fz2Values previousValues, Fz2Values currentValues,
			double Pe, double Rp, double C0, double beta){
		this.previousValues = previousValues;
		this.currentValues = currentValues;
		this.Pe = Pe;
		this.Rp = Rp;
		this.C0 = C0;
		this.beta = beta;
	}

	@Override
	public void run() {
		init();
		AtomicDouble totalChange = new AtomicDouble();
		int iterations =0;
		do{
			currentValues.p = new double[stepsX][stepsZ];
			totalChange.set(0);
			
			iterateXExcludeBordersParallel((x,jj) -> {
				iterateZExcludeBorders((z,m) -> {
					int j = jj;
					double[][] p = previousValues.p;
					double[][] c = previousValues.c;
					double[][] k = currentValues.k;
					double kk = k[j][m];
					double cc = c[j][m];
					double dkx = firstDerX(k, j, m);
					double dkz = firstDerZ(k, j, m);
					double dcx = firstDerX(c, j, m);
					double d2cx = secondDerX(c, j, m);
					double dcz = firstDerZ(c, j, m);
					double d2cz = secondDerZ(c, j, m);
					
					
					double a1 =  dkx/kk;
					double a2 =  - Rp/Pe*(1/kk*dkx*dcx + d2cx + 1/kk*dkz*dcz + d2cz);
					double a3 = + (p[j+1][m]+p[j-1][m])/dx/dx;
					double a4 = + (p[j][m+1]+p[j][m-1])/dz/dz;
					double val = a1 + a2 + a3 + a4;
					
					// REPLACED BY MORE SIMPLE CASE
//					double a1 = 1/kk*dkx;
//					double a2 = Rp/Pe*( 1/kk*z*dkx*dcx + z*d2cx + cc/kk*dkz + 2*dcz +z/kk*dkz*dcz + z*d2cz);
//					//a2 = a2 - Rp/Pe*dkz/kk/beta/C0; // EXCLUDED FROM EQUATIONS
//					double a3 = - dkx/kk * firstDerX(p, j, m) - dkz/kk * firstDerZ(p, j, m);
//					double a4 = - (p[j+1][m]+p[j-1][m])/dx/dx;
//					double a5 = - (p[j][m+1]+p[j][m-1])/dz/dz;
//					double val = a1 + a2 + a3 + a4 + a5;
					
					double cur = currentValues.p[j][m] = dz*dz*dx*dx/2/(dz*dz+dx*dx) * val;
					double prev = p[j][m];
					totalChange.addAndGet(Math.abs(cur - prev));
				});
			});
			iterateX((x,j) -> {
				if(j == 0)
					leftBoundaryCondition();
				else if(j == lastXInd)
					rightBoundaryCondition();
				else {
					bottomBoundaryCondition(j);
					topBoundaryCondition(j);
				}
			});
			
			previousValues.p = currentValues.p;
			iterations++;
		}
		while(totalChange.get() > stopCriteria && (iterations < 1_000 || totalChange.get() > stopCriteria*5));
		System.out.println(iterations+" "+totalChange.get());
	}

	private void leftBoundaryCondition() {
		currentValues.p[0] = new double[stepsZ];
//		iterateZ((z,m) -> {
//			currentValues.p[0][m] = currentValues.p[1][m] ; 	
//		});
		//Arrays.fill(currentValues.p[0], 1);
	}
	
	private void rightBoundaryCondition() {
		currentValues.p[lastXInd] = new double[stepsZ];
//		iterateZ((z,m) -> {
//			currentValues.p[lastXInd][m] = currentValues.p[lastXInd-1][m] ; 	
//		});
	}
	
	private void bottomBoundaryCondition(int j) {
		double[][] c = previousValues.c;
		currentValues.p[j][0] = currentValues.p[j][1] - dz*Rp/Pe*(c[j][0]);
		if(!simpleCase){
			currentValues.p[j][0] += -dz*Rp/Pe*0*firstDerZ(c, j, lastZInd);
		}
	}
	
	private void topBoundaryCondition(int j) {
		double[][] c = previousValues.c;
		currentValues.p[j][lastZInd] = dz*Rp/Pe*(c[j][lastZInd] + 1) - currentValues.p[j][lastZInd-1];
		if(!simpleCase){
			currentValues.p[j][lastZInd] += dz*Rp/Pe*maxZ*firstDerZ(c, j, lastZInd);
		}
	}
}
