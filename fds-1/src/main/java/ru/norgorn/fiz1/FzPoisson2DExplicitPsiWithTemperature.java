package ru.norgorn.fiz1;

public class FzPoisson2DExplicitPsiWithTemperature extends FzPoisson2DExplicit {
	
	private final double Ra;
	private final double Le;

	public FzPoisson2DExplicitPsiWithTemperature(Fz2Values previousValues, Fz2Values currentValues, double Pe, double Rp,
			double Ra, double Le, double C0) {
		super(previousValues, currentValues, Pe, Rp, C0);
		this.Ra = Ra;
		this.Le = Le;
	}
	
	@Override
	protected double solve(Integer j, Integer m) {
		double[][] p = previousValues.p;
		double[][] c = previousValues.c;
		double[][] k = currentValues.k;
		double kk = k[j][m];
		double dcx = firstDerX(c, j, m);
		double dcz = firstDerZ(c, j, m);
		double dtx = firstDerX(c, j, m);
		double dtz = firstDerZ(c, j, m);
		double dkx = firstDerX(k, j, m);
		double dkz = firstDerZ(k, j, m);
		
		double px = (p[j+1][m]+p[j-1][m])/dx;
		double pz = (p[j][m+1]+p[j][m-1])/dz;
		double dpx = firstDerX(p, j, m);
		double dpz = firstDerZ(p, j, m);
		
		double a1 = - Rp*kk*(dcz - dcx);
		double a2 = + px/dx;
		double a3 = + pz/dz;
		double a4 = (dkx*dpx + dkz*dpz)/kk;
		double a5 = - Le*Ra*kk*(dtz - dtx);
		double val = a1 + a2 + a3 + a4 + a5;
		
		double cur = dz*dz*dx*dx/2/(dz*dz+dx*dx) * val;
		return cur;
	}

	protected void leftBoundaryCondition() {
		currentValues.p[0] = new double[stepsZ];
		iterateZ((z, m) -> {
			currentValues.p[0][m] = Pe*z;
		});
	}
	
	protected void rightBoundaryCondition() {
		currentValues.p[lastXInd] = new double[stepsZ];
	}
	
	protected void bottomBoundaryCondition(int j) {
		if(j ==0)
			currentValues.p[j][0] = 0;
		else
			currentValues.p[j][0] = currentValues.p[j-1][0];
	}
	
	protected void topBoundaryCondition(int j) {
		if(j==0)
			currentValues.p[j][lastZInd] = Pe;
		else
			currentValues.p[j][lastZInd] = currentValues.p[j-1][lastZInd];
	}
}