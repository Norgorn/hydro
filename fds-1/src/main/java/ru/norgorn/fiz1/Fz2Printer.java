package ru.norgorn.fiz1;

import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.stream.Collectors;

public class Fz2Printer extends FzIterative {
	
	public static DecimalFormat formatGeneral = new DecimalFormat("#0.000");
	public static DecimalFormat formatPrecise = new DecimalFormat("#0.00000000");
	
	public static String formatRow(double[] c, DecimalFormat format) {
		String row = Arrays.stream(c).boxed()
				.map(d ->format.format(d)).collect(Collectors.joining("\t"))
				.replaceAll(",", "\\.");
		return row;
	}
	
	public static Fz2Printer clearFileAndPrepareWriter(String basePath, String parametersString, String name){
		return new Fz2Printer(basePath, parametersString, name);
	}
	
	private final String basePath;
	private final String parametersString;
	private final PrintWriter output;
	
	private volatile long startTime;	
	
	private Fz2Printer(String basePath, String parametersString, String name){
		try{
			this.basePath = basePath;
			this.parametersString = parametersString;
			String path = savePath(name);
			new PrintWriter(path).close();
			output = new PrintWriter(path);
			startTime = System.nanoTime();
			System.out.println(path);
		} catch(Exception e){
			throw new RuntimeException(e);
		}
	}
	
	public void printRow(double[] c, int i) {
		printRow(c, i, formatGeneral);
	}
	
	public void printRow(double[] c, int i, DecimalFormat format) {
		String row = formatRow(c, format);
		long end = System.nanoTime();
		long time = (end - startTime)/1_000_000;
		System.out.println(i+"("+time+"): "+row);
		startTime = end;
	}
	
	public void saveRow(double[][] row) {
		saveRow(row, formatGeneral);
	}
	
	public void saveVectorRow(double[][] vx, double[][] vz) {
		String[][] matrix = new String[stepsX][stepsZ];
		iterateX((x,j)-> {
			iterateZ((z,m)-> {
				matrix[j][m] = "{"+formatGeneral.format(vx[j][m]).replaceAll(",", "\\.")
						+","+formatGeneral.format(vz[j][m]).replaceAll(",", "\\.")+"}"; 
			});
		});
		
		String rowStr = Arrays.stream(matrix)
			.map(r ->Arrays.stream(r).collect(Collectors.joining(",")))
			.collect(Collectors.joining(","));
		
		print(rowStr);
	}
	
	public void saveRow(double[][] row, DecimalFormat format) {
		StringBuilder builder = new StringBuilder();
		for (int i=0; i< row.length; i++){
			builder.append(formatRow(row[i], format)).append(" ");
		}
		print(builder.toString());
	}
	
	private String savePath(String name) {
		return basePath+name+"_"+parametersString+".txt";
	}

	private void print(String str) {
		try {
			output.println(str);
			output.flush();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
