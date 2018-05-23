package ru.norgorn.fiz1;

import java.util.Arrays;

import com.google.common.util.concurrent.AtomicDouble;

public class FzPoisson2DExplicit extends FzDerivatives implements Runnable  {
	
	protected double stopCriteria = 0.0001;
	
	public final Fz2Values currentValues;
	protected final Fz2Values previousValues;
	protected final double Pe;
	protected final double Rp;
	protected final double C0;
	protected final double L=1;
	
	private final boolean simpleCase = true;
	
	public FzPoisson2DExplicit(Fz2Values previousValues, Fz2Values currentValues,
			double Pe, double Rp, double C0){
		this.previousValues = previousValues;
		this.currentValues = currentValues;
		this.Pe = Pe;
		this.Rp = Rp;
		this.C0 = C0;
	}

	@Override
	public void run() {
		init();
		AtomicDouble totalChange = new AtomicDouble();
		int iterations =0;
		do{
			currentValues.p = new double[stepsX][stepsZ];
			totalChange.set(0);
			
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
			iterateXExcludeBordersParallel((x,j) -> {
				iterateZExcludeBorders((z,m) -> {
					double prev = previousValues.p[j][m];
					double cur = solve(j, m);
					currentValues.p[j][m] = cur;
					totalChange.addAndGet(Math.abs(cur - prev));
				});
			});
			
			// Right boundary condition
			currentValues.p[lastXInd] = currentValues.p[lastXInd - 1];
			
			previousValues.p = currentValues.p;
			iterations++;
		}
		while(totalChange.get() > stopCriteria && (iterations < 5_000 || totalChange.get() > stopCriteria*5));
		System.out.println(iterations+" "+totalChange.get());
	}

	protected double solve(Integer j, Integer m) {
		double[][] p = previousValues.p;
		double[][] c = previousValues.c;
		double[][] k = currentValues.k;
		double kk = k[j][m];
		double cc = c[j][m];
		double dkx = firstDerX(k, j, m);
		double dkz = firstDerZ(k, j, m);
		double dcx = firstDerX(c, j, m);
		double dcz = firstDerZ(c, j, m);
		
		
		double a1 = 0;
		double a2 =  - Rp/Pe*(1/kk*dkx*cc + dcx + 1/kk*dkz*cc + dcz);
		double a3 = + (p[j+1][m]+p[j-1][m])/dx/dx;
		double a4 = + (p[j][m+1]+p[j][m-1])/dz/dz;
		double val = a1 + a2 + a3 + a4;
		
		double cur = dz*dz*dx*dx/2/(dz*dz+dx*dx) * val;
		return cur;
	}

	protected void leftBoundaryCondition() {
		currentValues.p[0] = new double[stepsZ];
		Arrays.fill(currentValues.p[0], 1);
	}
	
	protected void rightBoundaryCondition() {
		currentValues.p[lastXInd] = new double[stepsZ];
	}
	
	protected void bottomBoundaryCondition(int j) {
		double[][] c = previousValues.c;
		currentValues.p[j][0] = currentValues.p[j][1] - dz*Rp/Pe*(c[j][0]);
		if(!simpleCase){
			currentValues.p[j][0] += -dz*Rp/Pe*0*firstDerZ(c, j, lastZInd);
		}
	}
	
	protected void topBoundaryCondition(int j) {
		double[][] c = previousValues.c;
		currentValues.p[j][lastZInd] = dz*Rp/Pe*(c[j][lastZInd] + 1) - currentValues.p[j][lastZInd-1];
		if(!simpleCase){
			currentValues.p[j][lastZInd] += dz*Rp/Pe*maxZ*firstDerZ(c, j, lastZInd);
		}
	}
}
