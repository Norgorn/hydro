package ru.norgorn.fiz1;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.BiConsumer;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class Fz2MainIntegrationalValuesCounter {

	public static void main(String[] args) {
		try {
			new Fz2MainIntegrationalValuesCounter().run();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void run() throws Exception{
		int xSteps = 50;
		int zSteps = 50;
		int timeSteps = 999;
		
		
		List<String> paramStrings = Arrays.asList(
//				"C0=0.2_Fi=3.0_Q=5.0_a=15.0_b=19.0_Pe=100.0_Rp=10.0_beta=0.1"
//				,"C0=0.2_Fi=3.0_Q=5.0_a=15.0_b=19.0_Pe=100.0_Rp=15.0_beta=0.1"
//				,"C0=0.2_Fi=3.0_Q=5.0_a=15.0_b=19.0_Pe=100.0_Rp=20.0_beta=0.1"
//				,"C0=0.2_Fi=3.0_Q=5.0_a=15.0_b=19.0_Pe=100.0_Rp=25.0_beta=0.1"
//				"C0=0.2_Fi=3.0_Q=5.0_a=15.0_b=19.0_Pe=100.0_Rp=30.0_beta=0.1"
//				,"C0=0.2_Fi=3.0_Q=5.0_a=15.0_b=19.0_Pe=100.0_Rp=40.0_beta=0.1"
//				,"C0=0.2_Fi=3.0_Q=5.0_a=15.0_b=19.0_Pe=100.0_Rp=45.0_beta=0.1"
//				"C0=0.2_Fi=3.0_Q=5.0_a=15.0_b=19.0_Pe=100.0_Rp=50.0_beta=0.1"
//				,"C0=0.2_Fi=3.0_Q=5.0_a=15.0_b=19.0_Pe=100.0_Rp=55.0_beta=0.1"
				
//				"C0=0.2_Fi=3.0_Q=5.0_a=15.0_b=19.0_Pe=10.0_Rp=30.0_beta=0.1"
//				"C0=0.2_Fi=3.0_Q=5.0_a=15.0_b=19.0_Pe=50.0_Rp=30.0_beta=0.1"
//				"C0=0.2_Fi=3.0_Q=5.0_a=15.0_b=19.0_Pe=100.0_Rp=30.0_beta=0.1"
				"C0=0.2_Fi=3.0_Q=5.0_a=15.0_b=19.0_Pe=150.0_Rp=30.0_beta=0.1"
//				"C0=0.2_Fi=3.0_Q=5.0_a=15.0_b=19.0_Pe=200.0_Rp=30.0_beta=0.1"
//				"C0=0.2_Fi=3.0_Q=5.0_a=15.0_b=19.0_Pe=275.0_Rp=30.0_beta=0.1"
		);
		
		List<Fz2integrationalValues> results = IntStream.range(0, timeSteps+1)
				.mapToObj(i -> new Fz2integrationalValues())
				.collect(Collectors.toList());
		
		
		paramStrings.forEach(params -> {
			try {
				String cPath = "C:\\Users\\Sunny\\Documents\\Wolfram Mathematica\\2d_nd\\c_"+params+".txt";
				String qPath = "C:\\Users\\Sunny\\Documents\\Wolfram Mathematica\\2d_nd\\q_"+params+".txt";
				String vxPath = "C:\\Users\\Sunny\\Documents\\Wolfram Mathematica\\2d_nd\\vx_"+params+".txt";
				
				processLines(cPath, (i,values) -> {results.get(i).avgC = lastColumnAverage(xSteps, zSteps, values);});
				processLines(cPath, (i,values) -> {results.get(i).middleDropC = middleColumnDrop(xSteps, zSteps, values);});
				processLines(cPath, (i,values) -> {results.get(i).lastDropC = lastColumnDrop(xSteps, zSteps, values);});
				processLines(qPath, (i,values) -> {results.get(i).totalQ = totalSum(values);});
				processLines(vxPath, (i,values) -> {results.get(i).avgV = lastColumnAverage(xSteps, zSteps, values);});
			} catch (Exception e) {
				e.printStackTrace();
			};
		});
		
		
		System.out.println("avgC \t dropC_Middle \t dropC_Last \t totalQ \t avgV");
		AtomicInteger i = new AtomicInteger(0);
		results.forEach(result ->{
			if(i.getAndIncrement()%10 !=0)
				return;
			String str = result.avgC
					+" \t "+ result.middleDropC
					+" \t "+ result.lastDropC
					+" \t "+  result.totalQ
					+" \t "+  result.avgV;
			System.out.println(str.replaceAll("\\.", "\\,"));
		});
	}
	
	protected double lastColumnAverage(int xSteps, int zSteps, List<String> values) {
		return values.stream()
			.skip((xSteps-1)*zSteps)
			.mapToDouble(Double::parseDouble)
			.average().orElse(0);
	}

	protected double middleColumnDrop(int xSteps, int zSteps, List<String> values) {
		List<Double> middleColumn = values.stream()
			.skip((xSteps/2)*zSteps)
			.limit(zSteps)
			.map(Double::parseDouble)
			.collect(Collectors.toList());
		double d = middleColumn.get(zSteps-2) - middleColumn.get(2);
		return d;
	}
	
	protected double lastColumnDrop(int xSteps, int zSteps, List<String> values) {
		List<Double> middleColumn = values.stream()
			.skip((xSteps-1)*zSteps)
			.map(Double::parseDouble)
			.collect(Collectors.toList());
		double d = middleColumn.get(zSteps-2) - middleColumn.get(2);
		return d;
	}

	protected double totalSum(List<String> values) {
		return values.stream().mapToDouble(Double::parseDouble).sum();
//		return StreamEx.ofSubLists(values, xSteps)
//			.mapToDouble(column -> column.stream().mapToDouble(Double::parseDouble).sum())
//			.sum();
	}
	
	private void processLines(String filePath, BiConsumer<Integer, List<String>> valuesConsumer) throws Exception{
		BufferedReader buffer = new BufferedReader(new FileReader(filePath));
		String line;
		int i=0;
		while((line = buffer.readLine()) != null){
			List<String> values = Arrays.asList(line.split("\\s+"));
			valuesConsumer.accept(++i, values);
		}
		buffer.close();
	}
}
