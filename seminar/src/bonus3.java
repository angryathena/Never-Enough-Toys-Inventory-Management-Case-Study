import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Scanner;

import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

public class bonus3 {

	private static int y18[][];
	private static int y19[][];
	private static double revenue = 0;
	private static double relevance = 0;
	private static double serviceGT = 0;
	private static double serviceROT = 0;

	public static void main(String[] args) throws FileNotFoundException {
		// 2018 allocation under capacity AND service constraint
		set2018y();
		// Computing 2019 performance and setting allocation under capacity AND service
		// constraint
		set2019y();
		// Printing 2019 performance
//		PrintStream ps1 = new PrintStream(new File("bonus3performance2019.txt"));
//		ps1.println(revenue + "	" + relevance + "	" + serviceGT + "	" + serviceROT);
//		ps1.close();
//		// Computing 2020 performance and setting allocation under capacity NO service
//		// constraint
//		run2020();
//		PrintStream ps2 = new PrintStream(new File("bonus3performance2020.txt"));
//		ps2.println(revenue + "	" + relevance + "	" + serviceGT + "	" + serviceROT);
//		ps2.close();

	}

	public static void set2018y() throws FileNotFoundException {
		// Chunk list 2018
		Scanner chunkList = new Scanner(new File("chunks.txt"));
		ArrayList<String> chunks = new ArrayList<>();
		while (chunkList.hasNextLine()) {
			chunks.add(chunkList.nextLine());
		}
		chunkList.close();

		// Size list
		ArrayList<String> sizes = new ArrayList<>();
		sizes.add("3XS");
		sizes.add("XXS");
		sizes.add("XS");
		sizes.add("S");
		sizes.add("M");
		sizes.add("L");
		sizes.add("XL");
		sizes.add("XXL");
		sizes.add("3XL");

		// Warehouse list
		ArrayList<String> warehouses = new ArrayList<>();
		warehouses.add("Warehouse A");
		warehouses.add("Warehouse B");

		int S = 9;
		int K = chunks.size();
		int T = 52;
		int W = 2;

		// Demand list 2018
		Scanner demand = new Scanner(new File("demands.txt"));
		double[][][] d = new double[S][K][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				Scanner colDemand = new Scanner(demand.nextLine());
				for (int t = 0; t < T; t++) {
					d[s][k][t] = colDemand.nextDouble();
				}
				colDemand.close();
			}
		}
		demand.close();

		// Revenue list 2018
		Scanner profit = new Scanner(new File("revenues.txt"));
		double[][][] p = new double[S][K][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				Scanner colProfit = new Scanner(profit.nextLine());
				for (int t = 0; t < T; t++) {
					p[s][k][t] = colProfit.nextDouble();
				}
				colProfit.close();
			}
		}
		profit.close();

		// Volume list 2018
		Scanner volume = new Scanner(new File("volumes.txt"));
		double[][][] v = new double[S][K][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				Scanner colVolume = new Scanner(volume.nextLine());
				for (int t = 0; t < T; t++) {
					v[s][k][t] = colVolume.nextDouble();
				}
				colVolume.close();
			}
		}
		volume.close();

		double[] Q = { 2250, 450 };

		// Solving 2018, y[K][W] stored
		try {

			solveModel18(chunks, S, K, T, W, p, d, v, Q);

		} catch (IloException e) {
			System.out.println("A Cplex exception occured: " + e.getMessage());
			e.printStackTrace();
		}
	}

	public static void set2019y() throws FileNotFoundException {

		// Chunk list 2018
		Scanner chunkList = new Scanner(new File("chunks.txt"));
		ArrayList<String> chunks = new ArrayList<>();
		while (chunkList.hasNextLine()) {
			chunks.add(chunkList.nextLine());
		}
		chunkList.close();

		// Chunk list 2019
		Scanner chunkList2019 = new Scanner(new File("chunks2019.txt"));
		ArrayList<String> chunks2019 = new ArrayList<>();
		while (chunkList2019.hasNextLine()) {
			chunks2019.add(chunkList2019.nextLine());
		}
		chunkList2019.close();
		int K19 = chunks2019.size();

		// Size list
		ArrayList<String> sizes = new ArrayList<>();
		sizes.add("3XS");
		sizes.add("XXS");
		sizes.add("XS");
		sizes.add("S");
		sizes.add("M");
		sizes.add("L");
		sizes.add("XL");
		sizes.add("XXL");
		sizes.add("3XL");

		// Warehouse list
		ArrayList<String> warehouses = new ArrayList<>();
		warehouses.add("Warehouse A");
		warehouses.add("Warehouse B");

		int S = 9;
		int T = 52;
		int W = 2;

		// Demand list 2019
		Scanner demand2019 = new Scanner(new File("demands2019a2.txt"));
		double[][][] d2019 = new double[S][K19][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K19; k++) {
				Scanner colDemand = new Scanner(demand2019.nextLine());
				for (int t = 0; t < T; t++) {
					d2019[s][k][t] = colDemand.nextDouble();
				}
				colDemand.close();
			}
		}
		demand2019.close();

		// Revenue list 2019
		Scanner rev2019 = new Scanner(new File("revenues2019a2.txt"));
		double[][][] p2019 = new double[S][K19][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K19; k++) {
				Scanner colProfit = new Scanner(rev2019.nextLine());
				for (int t = 0; t < T; t++) {
					p2019[s][k][t] = colProfit.nextDouble();
				}
				colProfit.close();
			}
		}
		rev2019.close();

		// Volume list 2019
		Scanner volume2019 = new Scanner(new File("volumes2019a2.txt"));
		double[][][] v2019 = new double[S][K19][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K19; k++) {
				Scanner colVolume = new Scanner(volume2019.nextLine());
				for (int t = 0; t < T; t++) {
					v2019[s][k][t] = colVolume.nextDouble();
				}
				colVolume.close();
			}
		}
		volume2019.close();

		// Relevance score list 2019
		Scanner relevance2019 = new Scanner(new File("relevances2019a2.txt"));
		double[][][] r2019 = new double[S][K19][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K19; k++) {
				Scanner colRel = new Scanner(relevance2019.nextLine());
				for (int t = 0; t < T; t++) {
					r2019[s][k][t] = colRel.nextDouble();
				}
				colRel.close();
			}
		}
		relevance2019.close();

		double[] Q = { 2250, 450 };

		try {

			solveModel19(chunks2019, chunks, S, K19, T, W, p2019, d2019, v2019, r2019, Q);
			// Printing 2019 y values
			PrintStream ps3 = new PrintStream(new File("yValues19.txt"));
			for (int k=0; k<K19; k++) {
				for (int w=0; w<W; w++) {
					ps3.print(y19[k][w] + " ");
				}
				ps3.println();
			}
			ps3.close();

		} catch (IloException e) {
			System.out.println("A Cplex exception occured: " + e.getMessage());
			e.printStackTrace();
		}
	}

	public static void run2020() throws FileNotFoundException {
		// Chunk list 2018
		Scanner chunkList = new Scanner(new File("chunks.txt"));
		ArrayList<String> chunks = new ArrayList<>();
		while (chunkList.hasNextLine()) {
			chunks.add(chunkList.nextLine());
		}
		chunkList.close();

		// Chunk list 2019
		Scanner chunkList2019 = new Scanner(new File("chunks2019.txt"));
		ArrayList<String> chunks2019 = new ArrayList<>();
		while (chunkList2019.hasNextLine()) {
			chunks2019.add(chunkList2019.nextLine());
		}
		chunkList2019.close();

		// Chunk list 2020
		Scanner chunkList2020 = new Scanner(new File("chunks2020.txt"));
		ArrayList<String> chunks2020 = new ArrayList<>();
		while (chunkList2020.hasNextLine()) {
			chunks2020.add(chunkList2020.nextLine());
		}
		chunkList2020.close();

		// Size list
		ArrayList<String> sizes = new ArrayList<>();
		sizes.add("3XS");
		sizes.add("XXS");
		sizes.add("XS");
		sizes.add("S");
		sizes.add("M");
		sizes.add("L");
		sizes.add("XL");
		sizes.add("XXL");
		sizes.add("3XL");

		// Warehouse list
		ArrayList<String> warehouses = new ArrayList<>();
		warehouses.add("Warehouse A");
		warehouses.add("Warehouse B");
		double[] Q = { 2250, 450 };

		int S = sizes.size();
		int K = chunks.size();
		int T = 52;
		int W = warehouses.size();

		// Demand list
		Scanner demand18 = new Scanner(new File("demand18.txt"));
		double[][][] d18 = new double[S][K][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				Scanner colDemand18 = new Scanner(demand18.nextLine());
				for (int t = 0; t < T; t++) {
					d18[s][k][t] = colDemand18.nextDouble();
				}
				colDemand18.close();
			}
		}
		demand18.close();
		Scanner demand19 = new Scanner(new File("demand19.txt"));
		double[][][] d19 = new double[S][K][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				Scanner colDemand19 = new Scanner(demand19.nextLine());
				for (int t = 0; t < T; t++) {
					d19[s][k][t] = colDemand19.nextDouble();
				}
				colDemand19.close();
			}
		}
		demand19.close();

		// Profit list
		Scanner profit18 = new Scanner(new File("revenues18.txt"));
		double[][][] p18 = new double[S][K][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				Scanner colProfit18 = new Scanner(profit18.nextLine());
				for (int t = 0; t < T; t++) {
					p18[s][k][t] = colProfit18.nextDouble();
				}
				colProfit18.close();
			}
		}
		profit18.close();
		Scanner profit19 = new Scanner(new File("revenues19.txt"));
		double[][][] p19 = new double[S][K][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				Scanner colProfit19 = new Scanner(profit19.nextLine());
				for (int t = 0; t < T; t++) {
					p19[s][k][t] = colProfit19.nextDouble();
				}
				colProfit19.close();
			}
		}
		profit19.close();

		// Volume list
		Scanner volume18 = new Scanner(new File("volumes18.txt"));
		double[][][] v18 = new double[S][K][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				Scanner colVolume18 = new Scanner(volume18.nextLine());
				for (int t = 0; t < T; t++) {
					v18[s][k][t] = colVolume18.nextDouble();
				}
				colVolume18.close();
			}
		}
		volume18.close();
		Scanner volume19 = new Scanner(new File("volumes19.txt"));
		double[][][] v19 = new double[S][K][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				Scanner colVolume19 = new Scanner(volume19.nextLine());
				for (int t = 0; t < T; t++) {
					v19[s][k][t] = colVolume19.nextDouble();
				}
				colVolume19.close();
			}
		}
		volume19.close();

		// Relevance score list
		Scanner relev18 = new Scanner(new File("relevance18.txt"));
		double[][][] r18 = new double[S][K][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				Scanner colRel18 = new Scanner(relev18.nextLine());
				for (int t = 0; t < T; t++) {
					r18[s][k][t] = colRel18.nextDouble();
				}
				colRel18.close();
			}
		}
		relev18.close();
		Scanner relev19 = new Scanner(new File("relevance19.txt"));
		double[][][] r19 = new double[S][K][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				Scanner colRel19 = new Scanner(relev19.nextLine());
				for (int t = 0; t < T; t++) {
					r19[s][k][t] = colRel19.nextDouble();
				}
				colRel19.close();
			}
		}
		relev19.close();

		// 2020 week 1 data

		// Demand list
		Scanner demand = new Scanner(new File("demand20.txt"));
		double[][][] d20 = new double[S][K][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				Scanner colDemand = new Scanner(demand.nextLine());
				d20[s][k][0] = colDemand.nextDouble();
				colDemand.close();
			}
		}
		demand.close();

		// Profit list
		Scanner profit = new Scanner(new File("revenues20.txt"));
		double[][][] p20 = new double[S][K][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				Scanner colProfit = new Scanner(profit.nextLine());
				p20[s][k][0] = colProfit.nextDouble();
				colProfit.close();
			}
		}
		profit.close();

		// Volume list
		Scanner volume = new Scanner(new File("volumes20.txt"));
		double[][][] v20 = new double[S][K][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				Scanner colVolume = new Scanner(volume.nextLine());
				v20[s][k][0] = colVolume.nextDouble();
				colVolume.close();
			}
		}
		volume.close();

		// Relevance scores
		double[][][] r20 = new double[S][K][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				if (r19[s][k][0] == 0)
					r20[s][k][0] = r18[s][k][0];
				else
					r20[s][k][0] = r19[s][k][0];
			}
		}

		// Forecasting rest of 2020
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				for (int t = 1; t < T; t++) {
					// demand 2020
					if (d18[s][k][t] == 0)
						d20[s][k][t] = d19[s][k][t];
					else if (d19[s][k][t] == 0)
						d20[s][k][t] = d18[s][k][t];
					else
						d20[s][k][t] = Math.round(d19[s][k][t] * d19[s][k][t] / d18[s][k][t]);
					// price 2020
					if (p19[s][k][t] == 0)
						p20[s][k][t] = p18[s][k][t];
					else
						p20[s][k][t] = p19[s][k][t];
					// volume 2020
					if (v19[s][k][t] == 0)
						v20[s][k][t] = v18[s][k][t];
					else
						v20[s][k][t] = v19[s][k][t];
					// relevance 2020
					if (r19[s][k][t] == 0)
						r20[s][k][t] = r18[s][k][t];
					else
						r20[s][k][t] = r19[s][k][t];

				}
			}
		}

		double[][] dAgg = new double[K][T];
		double[][] vAgg = new double[K][T];
		double[][] pAgg = new double[K][T];
		double[][] rAgg = new double[K][T];

		for (int k = 0; k < K; k++) {
			for (int t = 0; t < T; t++) {
				dAgg[k][t] = 0;
				pAgg[k][t] = 0;
				vAgg[k][t] = 0;
				for (int s = 0; s < S; s++) {
					dAgg[k][t] += d20[s][k][t];
					pAgg[k][t] += d20[s][k][t] * p20[s][k][t];
					vAgg[k][t] += d20[s][k][t] * v20[s][k][t];
					if (r20[s][k][t] != 0)
						rAgg[k][t] = r20[s][k][t];
				}
				if (dAgg[k][t] > 0) {
					pAgg[k][t] = pAgg[k][t] / dAgg[k][t];
					vAgg[k][t] = vAgg[k][t] / dAgg[k][t];
				}
			}
		}

		try {
			revenueModel(chunks2020, chunks2019, chunks, 1, K, T, W, pAgg, dAgg, vAgg, rAgg, Q);
		} catch (IloException e) {
			System.out.println("A Cplex exception occured: " + e.getMessage());
			e.printStackTrace();
		}
	}

	public static void solveModel18(ArrayList<String> chunks, int S, int K, int T, int W, double[][][] p,
			double[][][] d, double[][][] v, double[] Q) throws IloException {
		// Model created
		IloCplex cplex = new IloCplex();
		final double V_GT = 0.98;
		final double V_ROT = 0.95;
		cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1e-4);
		// Decision variables
		IloNumVar[][][] x = new IloNumVar[S][K][T]; // size
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				for (int t = 0; t < T; t++) {
					x[s][k][t] = cplex.intVar(0, (int) d[s][k][t]);
				}
			}
		}
		IloNumVar[][] y = new IloNumVar[K][W]; // size
		for (int k = 0; k < K; k++) {
			for (int w = 0; w < W; w++) {
				y[k][w] = cplex.boolVar();
			}
		}
		IloNumVar[][][] q = new IloNumVar[K][W][T]; // size
		for (int k = 0; k < K; k++) {
			for (int w = 0; w < W; w++) {
				for (int t = 0; t < T; t++) {
					q[k][w][t] = cplex.numVar(0, Q[w]);
				}
			}
		}

		// Objective function
		IloNumExpr objExpr = cplex.constant(0);
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				for (int t = 0; t < T; t++) {
					objExpr = cplex.sum(objExpr, cplex.prod(x[s][k][t], p[s][k][t]));
				}
			}
		}
		cplex.addMaximize(objExpr);

		// Single warehouse per chunk constraint
		IloNumExpr[] sum_y = new IloNumExpr[K];
		for (int k = 0; k < K; k++) {
			sum_y[k] = cplex.constant(0);
			for (int w = 0; w < W; w++) {
				sum_y[k] = cplex.sum(sum_y[k], y[k][w]);
			}
			// System.out.println(sum_y[k]);
			cplex.addEq(sum_y[k], 1);
		}
		// Capacity constraint
		for (int k = 0; k < K; k++) {
			for (int w = 0; w < W; w++) {
				for (int t = 0; t < T; t++) {
					cplex.addLe(q[k][w][t], cplex.prod(y[k][w], Q[w]));
				}
			}
		}

		IloNumExpr[][] sum_xv = new IloNumExpr[K][T];
		IloNumExpr[][] sum_w_q = new IloNumExpr[K][T];
		for (int t = 0; t < T; t++) {
			for (int k = 0; k < K; k++) {
				sum_xv[k][t] = cplex.constant(0);
				sum_w_q[k][t] = cplex.constant(0);
				for (int s = 0; s < S; s++) {
					sum_xv[k][t] = cplex.sum(sum_xv[k][t], cplex.prod(x[s][k][t], v[s][k][t]));
				}
				for (int w = 0; w < W; w++) {
					sum_w_q[k][t] = cplex.sum(sum_w_q[k][t], q[k][w][t]);
				}
				cplex.addEq(sum_xv[k][t], sum_w_q[k][t]);
			}
		}

		IloNumExpr[][] sum_k_q = new IloNumExpr[W][T];
		for (int w = 0; w < W; w++) {
			for (int t = 0; t < T; t++) {
				sum_k_q[w][t] = cplex.constant(0);
				for (int k = 0; k < K; k++) {
					sum_k_q[w][t] = cplex.sum(sum_k_q[w][t], q[k][w][t]);
				}
				cplex.addLe(sum_k_q[w][t], Q[w]);
			}
		}

		// Service Level Constraints
		IloNumExpr sum_xGT = cplex.constant(0);
		IloNumExpr sum_dGT = cplex.constant(0);
		IloNumExpr sum_xROT = cplex.constant(0);
		IloNumExpr sum_dROT = cplex.constant(0);
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				for (int t = 0; t < T; t++) {
					if (chunks.get(k).substring(0, 1).equals("R")) {
						sum_xROT = cplex.sum(sum_xROT, x[s][k][t]);
						sum_dROT = cplex.sum(sum_dROT, d[s][k][t]);

					} else {
						sum_xGT = cplex.sum(sum_xGT, x[s][k][t]);
						sum_dGT = cplex.sum(sum_dGT, d[s][k][t]);

					}

				}

			}

		}
		cplex.addGe(sum_xROT, cplex.prod(cplex.constant(V_ROT), sum_dROT));
		cplex.addGe(sum_xGT, cplex.prod(cplex.constant(V_GT), sum_dGT));

		// Engaging solver and formatting console output
		cplex.solve();
		if (cplex.getStatus() == IloCplex.Status.Optimal) {
			System.out.println("Found optimal solution!");
			System.out.println("Objective = " + cplex.getObjValue());
			y18 = new int[K][W];
			for (int k = 0; k < K; k++) {
				for (int w = 0; w < W; w++) {
					y18[k][w] = (int) cplex.getValue(y[k][w]);
				}
			}

		} else {
			System.out.println("No optimal solution found");
		}
		cplex.close();
	}

	public static void solveModel19(ArrayList<String> chunks, ArrayList<String> chunks18, int S, int K, int T, int W,
			double[][][] p, double[][][] d, double[][][] v, double[][][] r, double[] Q)
			throws IloException, FileNotFoundException {
		// Model created
		IloCplex cplex = new IloCplex();
		final double V_GT = 0.98;
		final double V_ROT = 0.95;
		cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1e-4);
		// Decision variables
		IloNumVar[][][] x = new IloNumVar[S][K][T]; // size
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				for (int t = 0; t < T; t++) {
					x[s][k][t] = cplex.intVar(0, (int) d[s][k][t]);
				}
			}
		}

		// Limited number of y variables needed
		IloNumVar[][] y = new IloNumVar[K][W]; // size
		for (int k = 0; k < K; k++) {
			for (int w = 0; w < W; w++) {
				y[k][w] = cplex.boolVar();
				if (chunks18.contains(chunks.get(k)))
					cplex.addEq(y[k][w], cplex.constant(y18[chunks18.indexOf(chunks.get(k))][w]));
			}
		}

		IloNumVar[][][] q = new IloNumVar[K][W][T]; // size
		for (int k = 0; k < K; k++) {
			for (int w = 0; w < W; w++) {
				for (int t = 0; t < T; t++) {
					q[k][w][t] = cplex.numVar(0, Q[w]);
				}
			}
		}

		// Objective function
		IloNumExpr objExpr = cplex.constant(0);
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				for (int t = 0; t < T; t++) {
					objExpr = cplex.sum(objExpr, cplex.prod(x[s][k][t], p[s][k][t]));
				}
			}
		}
		cplex.addMaximize(objExpr);

		// Single warehouse per chunk constraint
		IloNumExpr[] sum_y = new IloNumExpr[K];
		for (int k = 0; k < K; k++) {
			sum_y[k] = cplex.constant(0);
			for (int w = 0; w < W; w++) {
				sum_y[k] = cplex.sum(sum_y[k], y[k][w]);
			}
			if (!chunks18.contains(chunks.get(k)))
				cplex.addEq(sum_y[k], 1);
		}
		// Capacity constraint
		for (int k = 0; k < K; k++) {
			for (int w = 0; w < W; w++) {
				for (int t = 0; t < T; t++) {
					cplex.addLe(q[k][w][t], cplex.prod(y[k][w], Q[w]));
				}
			}
		}

		IloNumExpr[][] sum_xv = new IloNumExpr[K][T];
		IloNumExpr[][] sum_w_q = new IloNumExpr[K][T];
		for (int t = 0; t < T; t++) {
			for (int k = 0; k < K; k++) {
				sum_xv[k][t] = cplex.constant(0);
				sum_w_q[k][t] = cplex.constant(0);
				for (int s = 0; s < S; s++) {
					sum_xv[k][t] = cplex.sum(sum_xv[k][t], cplex.prod(x[s][k][t], v[s][k][t]));
				}
				for (int w = 0; w < W; w++) {
					sum_w_q[k][t] = cplex.sum(sum_w_q[k][t], q[k][w][t]);
				}
				cplex.addEq(sum_xv[k][t], sum_w_q[k][t]);
			}
		}

		IloNumExpr[][] sum_k_q = new IloNumExpr[W][T];
		for (int w = 0; w < W; w++) {
			for (int t = 0; t < T; t++) {
				sum_k_q[w][t] = cplex.constant(0);
				for (int k = 0; k < K; k++) {
					sum_k_q[w][t] = cplex.sum(sum_k_q[w][t], q[k][w][t]);
				}
				cplex.addLe(sum_k_q[w][t], Q[w]);
			}
		}

		// Service Level Constraints
		IloNumExpr sum_xGT = cplex.constant(0);
		IloNumExpr sum_dGT = cplex.constant(0);
		IloNumExpr sum_xROT = cplex.constant(0);
		IloNumExpr sum_dROT = cplex.constant(0);
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				for (int t = 0; t < T; t++) {
					if (chunks.get(k).substring(0, 1).equals("R")) {
						sum_xROT = cplex.sum(sum_xROT, x[s][k][t]);
						sum_dROT = cplex.sum(sum_dROT, d[s][k][t]);

					} else {
						sum_xGT = cplex.sum(sum_xGT, x[s][k][t]);
						sum_dGT = cplex.sum(sum_dGT, d[s][k][t]);

					}

				}

			}

		}
		cplex.addGe(sum_xROT, cplex.prod(cplex.constant(V_ROT), sum_dROT));
		cplex.addGe(sum_xGT, cplex.prod(cplex.constant(V_GT), sum_dGT));

		// Engaging solver and formatting console output
		cplex.solve();
		if (cplex.getStatus() == IloCplex.Status.Optimal) {
			System.out.println("Found optimal solution!");
			revenue = cplex.getObjValue();
			y19 = new int[K][W];
			for (int k = 0; k < K; k++) {
				for (int w = 0; w < W; w++) {
					y19[k][w] = (int) cplex.getValue(y[k][w]);
				}
			}
			double totalDemandGT = 0;
			double totalDemandROT = 0;
			double totalSalesGT = 0;
			double totalSalesROT = 0;
			for (int t = 0; t < T; t++) {
				for (int k = 0; k < K; k++) {
					for (int s = 0; s < S; s++) {
						if (chunks.get(k).substring(0, 1).equals("R")) {
							totalSalesROT += cplex.getValue(x[s][k][t]);
							totalDemandROT += d[s][k][t];

						} else {
							totalSalesGT += cplex.getValue(x[s][k][t]);
							totalDemandGT += d[s][k][t];
						}
						relevance += r[s][k][t] * (int) cplex.getValue(x[s][k][t]);

					}

				}

			}
			serviceGT = totalSalesGT / totalDemandGT;
			serviceROT = totalSalesROT / totalDemandROT;

		} else {
			System.out.println("No optimal solution found");
		}
		cplex.close();
	}

	public static void revenueModel(ArrayList<String> chunks, ArrayList<String> chunks19, ArrayList<String> chunks18,
			int S, int K, int T, int W, double[][] pAgg, double[][] dAgg, double[][] vAgg, double[][] rAgg, double[] Q)
			throws IloException {
		// Model created
		IloCplex cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1e-3);
		// cplex.setOut(null);

		// Decision variables
		IloNumVar[][] x = new IloNumVar[K][T]; // size
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				for (int t = 0; t < T; t++) {
					x[k][t] = cplex.intVar(0, Integer.MAX_VALUE);
				}
			}
		}
		IloNumVar[][] y = new IloNumVar[K][W]; // size
		for (int k = 0; k < K; k++) {
			for (int w = 0; w < W; w++) {
				y[k][w] = cplex.boolVar();
				if (chunks19.contains(chunks.get(k)))
					cplex.addEq(y[k][w], y19[chunks19.indexOf(chunks.get(k))][w]);
				else if (chunks18.contains(chunks.get(k)))
					cplex.addEq(y[k][w], y18[chunks18.indexOf(chunks.get(k))][w]);
			}
		}
		IloNumVar[][][] q = new IloNumVar[K][W][T]; // size
		for (int k = 0; k < K; k++) {
			for (int w = 0; w < W; w++) {
				for (int t = 0; t < T; t++) {
					q[k][w][t] = cplex.numVar(0, Q[w]);
				}
			}
		}
		IloNumVar[][] g = new IloNumVar[K][9]; // size

		for (int k = 0; k < K; k++) {
			for (int t = 0; t < 9; t++) {
				g[k][t] = cplex.boolVar();
			}
		}

		IloNumVar[][] z = new IloNumVar[K][9]; // size
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				for (int t = 0; t < 9; t++) {
					z[k][t] = cplex.intVar(0, (int) dAgg[k][t]);
				}
			}
		}
		IloNumVar[][] I = new IloNumVar[K][9];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				for (int t = 0; t < 9; t++) {
					I[k][t] = cplex.intVar(0, Integer.MAX_VALUE);
				}
			}
		}
		// Objective function: 4a
		IloNumExpr objExpr = cplex.constant(0);
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				for (int t = 0; t < 43; t++) {
					objExpr = cplex.sum(objExpr, cplex.prod(x[k][t], pAgg[k][t]));
					cplex.addLe(x[k][t], dAgg[k][t]);
				}
				for (int t = 43; t < T; t++) {
					objExpr = cplex.sum(objExpr, cplex.prod(z[k][t - 43], pAgg[k][t]));
				}
			}
		}
		cplex.addMaximize(objExpr);

		// Upper bounds on order size: 4b, 4c
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				for (int t = 43; t < 51; t++) {
					cplex.addLe(x[k][t], cplex.prod(dAgg[k][t] + dAgg[k][t + 1], g[k][t - 43]));
				}
				cplex.addLe(x[k][51], cplex.prod(dAgg[k][51], g[k][8]));
			}
		}

		// Order buffer: 4d
		for (int k = 0; k < K; k++) {
			for (int t = 0; t < 8; t++) {
				cplex.addLe(cplex.sum(g[k][t], g[k][t + 1]), 1);
			}
		}

		// Sales, orders and inventory: 4e

		for (int k = 0; k < K; k++) {
			cplex.addEq(I[k][0], cplex.diff(x[k][43], z[k][0]));
			for (int t = 44; t < T; t++) {
				cplex.addEq(I[k][t - 43], cplex.sum(I[k][t - 44], cplex.diff(x[k][t], z[k][t - 43])));
			}

		}

		// Single warehouse per chunk: 4f
		IloNumExpr[] sum_y = new IloNumExpr[K];
		for (int k = 0; k < K; k++) {
			sum_y[k] = cplex.constant(0);
			for (int w = 0; w < W; w++) {
				sum_y[k] = cplex.sum(sum_y[k], y[k][w]);
			}
			cplex.addEq(sum_y[k], 1);
		}
		// Capacity constraint: 4g, 4h, 4i
		for (int k = 0; k < K; k++) {
			for (int w = 0; w < W; w++) {
				for (int t = 0; t < T; t++) {
					cplex.addLe(q[k][w][t], cplex.prod(y[k][w], Q[w]));
				}
			}
		}
		IloNumExpr[][] sum_xv = new IloNumExpr[K][T];
		IloNumExpr[][] sum_w_q = new IloNumExpr[K][T];
		for (int t = 0; t < T; t++) {
			for (int k = 0; k < K; k++) {
				sum_xv[k][t] = cplex.constant(0);
				sum_w_q[k][t] = cplex.constant(0);
				if (t < 43) {
					sum_xv[k][t] = cplex.sum(sum_xv[k][t], cplex.prod(x[k][t], vAgg[k][t]));
				} else {
					sum_xv[k][t] = cplex.sum(sum_xv[k][t], cplex.prod(I[k][t - 43], vAgg[k][t]));
				}
				for (int w = 0; w < W; w++) {
					sum_w_q[k][t] = cplex.sum(sum_w_q[k][t], q[k][w][t]);
				}
				cplex.addEq(sum_xv[k][t], sum_w_q[k][t]);
			}
		}
		IloNumExpr[][] sum_k_q = new IloNumExpr[W][T];
		for (int w = 0; w < W; w++) {
			for (int t = 0; t < T; t++) {
				sum_k_q[w][t] = cplex.constant(0);
				for (int k = 0; k < K; k++) {
					sum_k_q[w][t] = cplex.sum(sum_k_q[w][t], q[k][w][t]);
				}
				cplex.addLe(sum_k_q[w][t], Q[w]);
			}
		}

		// Engaging solver and formatting console output
		cplex.solve();
		if (cplex.getStatus() == IloCplex.Status.Optimal) {
			revenue = cplex.getObjValue();
			relevance = 0;
			double totalDemandGT = 0;
			double totalDemandROT = 0;
			double totalSalesGT = 0;
			double totalSalesROT = 0;
			for (int k = 0; k < K; k++) {
				for (int t = 0; t < T; t++) {
					relevance += rAgg[k][t] * (int) cplex.getValue(x[k][t]);
					if (chunks.get(k).substring(0, 1).equals("R")) {
						totalSalesROT += cplex.getValue(x[k][t]);
						totalDemandROT += dAgg[k][t];

					} else {
						totalSalesGT += cplex.getValue(x[k][t]);
						totalDemandGT += dAgg[k][t];
					}
				}
			}

			serviceGT = totalSalesGT / totalDemandGT;
			serviceROT = totalSalesROT / totalDemandROT;

		} else
			System.out.println("No optimal solution found in maxRevenue");
		cplex.close();
	}
}
