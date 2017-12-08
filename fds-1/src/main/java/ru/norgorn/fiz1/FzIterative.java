package ru.norgorn.fiz1;

import java.util.function.BiConsumer;
import java.util.stream.IntStream;

public class FzIterative {
	public int stepsT = 200_000;
	public int stepsX = 50;
	public int stepsZ = 50;
	public double maxT = 0.1;
	public double maxX = 1;
	public double maxZ = 1;
	
	protected int lastTInd;
	protected int lastXInd;
	protected int lastZInd;
	protected double dt;
	protected double dx;
	protected double dz;
	
	public void init(){
		lastTInd = lastTInd();
		lastXInd = lastXInd();
		lastZInd = lastZInd();
		dt = dt();
		dx = dx();
		dz = dz();
	}
	
	protected int lastTInd() {
		return stepsT-1;
	}
	
	protected int lastXInd() {
		return stepsX-1;
	}
	
	protected int lastZInd() {
		return stepsZ-1;
	}
	
	protected double dt() {
		return maxT/stepsT;
	}
	
	protected double dx() {
		return maxX/stepsX;
	}
	
	protected double dz() {
		return maxZ/stepsZ;
	}
	
	public void iterateT(BiConsumer<Double, Integer> operation){
		iterate(dt(), lastTInd(), operation);
	}
	
	public void iterateTParallel(BiConsumer<Double, Integer> operation){
		iterateParallel(dt(), lastTInd(), operation);
	}
	
	public void iterateX(BiConsumer<Double, Integer> operation){
		iterate(dx(), lastXInd(), operation);
	}
	
	public void iterateXParallel(BiConsumer<Double, Integer> operation){
		iterateParallel(dx(), lastXInd(), operation);
	}
	
	public void iterateXExcludeBorders(BiConsumer<Double, Integer> operation){
		iterateExcludeBorders(dx(), lastXInd(), operation);
	}
	
	public void iterateBackwardXExcludeBorders(BiConsumer<Double, Integer> operation){
		iterateBackwardExcludeBorders(dx(), lastXInd(), operation);
	}
	
	public void iterateXExcludeBordersParallel(BiConsumer<Double, Integer> operation){
		iterateExcludeBordersParallel(dx(), lastXInd(), operation);
	}
	
	public void iterateZ(BiConsumer<Double, Integer> operation){
		iterate(dz(), lastZInd(), operation);
	}
	
	public void iterateBackwardZ(BiConsumer<Double, Integer> operation){
		iterateBackward(dz(), lastZInd(), operation);
	}
	
	public void iterateZParallel(BiConsumer<Double, Integer> operation){
		iterateParallel(dz(), lastZInd(), operation);
	}
	
	public void iterateZExcludeBorders(BiConsumer<Double, Integer> operation){
		iterateExcludeBorders(dz(), lastZInd(), operation);
	}
	
	public void iterateBackwardZExcludeBorders(BiConsumer<Double, Integer> operation){
		iterateBackwardExcludeBorders(dz(), lastZInd(), operation);
	}
	
	public void iterateZExcludeBordersParallel(BiConsumer<Double, Integer> operation){
		iterateExcludeBordersParallel(dz(), lastZInd(), operation);
	}
	
	private static void iterate(double step, int lastIndex, BiConsumer<Double, Integer> operation){
		double x=0;
		for(int i=0; i<= lastIndex; i++){
			operation.accept(x, i);
			x+=step;
		}
	}
	
	private static void iterateParallel(double step, int lastIndex, BiConsumer<Double, Integer> operation){
		IntStream.range(0, lastIndex+1).parallel()
		.forEach(m -> {
			double z = m*step;
			operation.accept(z, m);
		});
	}
	
	private static void iterateExcludeBorders(double step, int lastIndex, BiConsumer<Double, Integer> operation){
		double x=step;
		for(int i=1; i< lastIndex; i++){
			operation.accept(x, i);
			x+=step;
		}
	}
	
	private static void iterateBackward(double step, int lastIndex, BiConsumer<Double, Integer> operation){
		double x=step*lastIndex;
		for(int i=lastIndex; i>= 0; i--){
			operation.accept(x, i);
			x+=step;
		}
	}
	
	private static void iterateBackwardExcludeBorders(double step, int lastIndex, BiConsumer<Double, Integer> operation){
		double x=step*(lastIndex-1);
		for(int i=lastIndex-1; i> 0; i--){
			operation.accept(x, i);
			x+=step;
		}
	}
	
	private static void iterateExcludeBordersParallel(double step, int lastIndex, BiConsumer<Double, Integer> operation){
		IntStream.range(1, lastIndex).parallel()
		.forEach(m -> {
			double z = m*step;
			operation.accept(z, m);
		});
	}
}
