import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Scanner;

import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

public class bonus2 {

	public static double maxRevenue = 1;
	public static double maxRelevance = 1;
	public static double maxRevenue2019 = 1;
	public static double maxRelevance2019 = 1;
	private static double serviceGT;
	public static double serviceROT;
	public static double revenue = 0;
	public static double relevance = 0;
	public static double objective = 0;
	public static double totalSalesROT = 0;
	public static double totalDemandROT = 0;
	public static double totalSalesGT = 0;
	public static double totalDemandGT = 0;
	public static int[][] yEq = new int[495][2];
	public static int[][][] xEq = new int[9][495][52];
	public static double[][][] qEq = new double[495][2][52];

	public static void main(String[] args) throws FileNotFoundException {
		// Chunk list
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

		// Demand list
		Scanner demand = new Scanner(new File("demands.txt"));
		double[][][] d = new double[S][K][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				Scanner colDemand = new Scanner(demand.nextLine());
				for (int t = 0; t < T; t++) {
					d[s][k][t] = colDemand.nextDouble();
					if (t < 45 || t > 47)
						d[s][k][t] = 0;
				}
				colDemand.close();
			}
		}
		demand.close();

		// Revenue list
		Scanner profit = new Scanner(new File("revenues.txt"));
		double[][][] p = new double[S][K][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				Scanner colProfit = new Scanner(profit.nextLine());
				for (int t = 0; t < T; t++) {
					p[s][k][t] = colProfit.nextDouble();
					if (t < 45 || t > 47)
						p[s][k][t] = 0;
				}
				colProfit.close();
			}
		}
		profit.close();

		// Volume list
		Scanner volume = new Scanner(new File("volumes.txt"));
		double[][][] v = new double[S][K][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				Scanner colVolume = new Scanner(volume.nextLine());
				for (int t = 0; t < T; t++) {
					v[s][k][t] = colVolume.nextDouble();
					if (t < 45 || t > 47)
						v[s][k][t] = 0;
				}
				colVolume.close();
			}
		}
		volume.close();

		// Relevance score list
		Scanner relev = new Scanner(new File("relevances.txt"));
		double[][][] r = new double[S][K][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				Scanner colRel = new Scanner(relev.nextLine());
				for (int t = 0; t < T; t++) {
					r[s][k][t] = colRel.nextDouble();
					if (t < 45 || t > 47)
						r[s][k][t] = 0;
				}
				colRel.close();
			}
		}
		relev.close();

		// Chunk list 2019
		Scanner chunkList2019 = new Scanner(new File("chunks2019.txt"));
		ArrayList<String> chunks2019 = new ArrayList<>();
		while (chunkList2019.hasNextLine()) {
			chunks2019.add(chunkList2019.nextLine());
		}
		chunkList2019.close();
		int K19 = chunks2019.size();

		// Demand list 2019
		Scanner demand2019 = new Scanner(new File("demands2019a2.txt"));
		double[][][] d2019 = new double[S][K19][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K19; k++) {
				Scanner colDemand = new Scanner(demand2019.nextLine());
				for (int t = 0; t < T; t++) {
					d2019[s][k][t] = colDemand.nextDouble();
					if (t < 45 || t > 47)
						d2019[s][k][t] = 0;
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
					if (t < 45 || t > 47)
						p2019[s][k][t] = 0;
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
					if (t < 45 || t > 47)
						v2019[s][k][t] = 0;
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
					if (t < 45 || t > 47)
						r2019[s][k][t] = 0;
				}
				colRel.close();
			}
		}
		relevance2019.close();

		double[] Q = { 2250, 450 };

		try {

			revenueModel(chunks2019, S, K19, T, W, p2019, d2019, v2019, r2019, Q);
			relevanceModel(chunks2019, S, K19, T, W, p2019, d2019, v2019, r2019, Q);
			maxRevenue2019 = maxRevenue;
			maxRelevance2019 = maxRelevance;
//						 maxRevenue2019 = 5075616.4027;
//						 maxRelevance2019 = 153493.069;
			revenueModel(chunks, S, K, T, W, p, d, v, r, Q);
			relevanceModel(chunks, S, K, T, W, p, d, v, r, Q);
//						 maxRevenue = 3646406.53;
//						 maxRelevance = 106048.245;

			PrintStream ps5 = new PrintStream(new File("bonus2AlphaResults.txt"));
			for (double alpha = 0; alpha <= 1; alpha += 0.1) {

				totalSalesROT = 0;
				totalDemandROT = 0;
				totalSalesGT = 0;
				totalDemandGT = 0;

				solveModel(chunks, S, K, T, W, p, d, r, v, Q, alpha, maxRevenue, maxRelevance);
				ps5.println(
						alpha + "	"+ revenue + "	" + relevance + "	" + totalSalesGT/totalDemandGT + "	" + totalSalesROT/totalDemandROT);
				
			}
			ps5.println("2019:");
			for (double alpha = 0; alpha <= 1; alpha += 0.1) {

				totalSalesROT = 0;
				totalDemandROT = 0;
				totalSalesGT = 0;
				totalDemandGT = 0;

				solveModel(chunks2019, S, K19, T, W, p2019, d2019, r2019, v2019, Q, alpha, maxRevenue2019,
						maxRelevance2019);
				ps5.println(
						alpha + "	"+ revenue + "	" + relevance + "	" + totalSalesGT/totalDemandGT + "	" + totalSalesROT/totalDemandROT);

			}

			ps5.close();

			System.out.println("---------------------------------------------------------------------------");

		} catch (IloException e) {
			System.out.println("A Cplex exception occured: " + e.getMessage());
			e.printStackTrace();
		}
	}


	public static void solveModel(ArrayList<String> chunks, int S, int K, int T, int W, double[][][] p, double[][][] d,
			double[][][] r, double[][][] v, double[] Q, double alpha, double maxRevenue, double maxRelevance)
			throws IloException {
		// Model created
		IloCplex cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1e-6);
		final double V_GT = 0.98;
		final double V_ROT = 0.95;
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
		IloNumExpr objExprRevenue = cplex.constant(0);
		IloNumExpr objExprRelevance = cplex.constant(0);
		IloNumExpr objFinal = cplex.constant(0);
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				for (int t = 0; t < T; t++) {
					objExprRevenue = cplex.sum(objExprRevenue, cplex.prod(x[s][k][t], p[s][k][t]));
					objExprRelevance = cplex.sum(objExprRelevance, cplex.prod(x[s][k][t], r[s][k][t]));
				}
			}
		}
		objFinal = cplex.sum(cplex.prod(cplex.constant(alpha / maxRevenue), objExprRevenue),
				cplex.prod(cplex.constant((1 - alpha) / maxRelevance), objExprRelevance));
		cplex.addMaximize(objFinal);

		// Service Level Constraints
//		IloNumExpr sum_xGT = cplex.constant(0);
//		IloNumExpr sum_dGT = cplex.constant(0);
//		IloNumExpr sum_xROT = cplex.constant(0);
//		IloNumExpr sum_dROT = cplex.constant(0);
//		for (int s = 0; s < S; s++) {
//			for (int k = 0; k < K; k++) {
//				for (int t = 0; t < T; t++) {
//					if (chunks.get(k).substring(0, 1).equals("R")) {
//						sum_xROT = cplex.sum(sum_xROT, x[s][k][t]);
//						sum_dROT = cplex.sum(sum_dROT, d[s][k][t]);
//
//					} else {
//						sum_xGT = cplex.sum(sum_xGT, x[s][k][t]);
//						sum_dGT = cplex.sum(sum_dGT, d[s][k][t]);
//
//					}
//
//				}
//
//			}
//
//		}
//		cplex.addGe(sum_xROT, cplex.prod(cplex.constant(V_ROT), sum_dROT));
//		cplex.addGe(sum_xGT, cplex.prod(cplex.constant(V_GT), sum_dGT));

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

		// Engaging solver and formatting console output
		cplex.solve();
		if (cplex.getStatus() == IloCplex.Status.Optimal) {
			System.out.println("Found optimal solution!");
			System.out.println("Objective = " + cplex.getObjValue());
			objective = cplex.getObjValue();
			revenue = cplex.getValue(objExprRevenue);
			relevance = cplex.getValue(objExprRelevance);

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
						xEq[s][k][t] = (int) cplex.getValue(x[s][k][t]);
					}
				}
			}

			for (int k = 0; k < K; k++) {
				for (int w = 0; w < W; w++) {
					yEq[k][w] = (int) cplex.getValue(y[k][w]);
				}
			}
			for (int k = 0; k < K; k++) {
				for (int w = 0; w < W; w++) {
					for (int t = 0; t < T; t++) {
						qEq[k][w][t] = cplex.getValue(q[k][w][t]);
					}

				}
			}
		} else {
			System.out.println("No optimal solution found");
		}
		cplex.close();
	}

	public static void revenueModel(ArrayList<String> chunks, int S, int K, int T, int W, double[][][] p,
			double[][][] d, double[][][] v, double[][][] r, double[] Q) throws IloException {
		// Model created
		final double V_GT = 0.98;
		final double V_ROT = 0.95;
		IloCplex cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1e-6);
		// cplex.setOut(null);

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
//		IloNumExpr sum_xGT = cplex.constant(0);
//		IloNumExpr sum_dGT = cplex.constant(0);
//		IloNumExpr sum_xROT = cplex.constant(0);
//		IloNumExpr sum_dROT = cplex.constant(0);
//		for (int s = 0; s < S; s++) {
//			for (int k = 0; k < K; k++) {
//				for (int t = 0; t < T; t++) {
//					if (chunks.get(k).substring(0, 1).equals("R")) {
//						sum_xROT = cplex.sum(sum_xROT, x[s][k][t]);
//						sum_dROT = cplex.sum(sum_dROT, d[s][k][t]);
//
//					} else {
//						sum_xGT = cplex.sum(sum_xGT, x[s][k][t]);
//						sum_dGT = cplex.sum(sum_dGT, d[s][k][t]);
//
//					}
//
//				}
//
//			}
//
//		}
//		cplex.addGe(sum_xROT, cplex.prod(cplex.constant(V_ROT), sum_dROT));
//		cplex.addGe(sum_xGT, cplex.prod(cplex.constant(V_GT), sum_dGT));

		// Engaging solver and formatting console output
		cplex.solve();
		if (cplex.getStatus() == IloCplex.Status.Optimal)
			maxRevenue = cplex.getObjValue();
		else
			System.out.println("No optimal solution found in maxRevenue");
		cplex.close();
	}

	public static void relevanceModel(ArrayList<String> chunks, int S, int K, int T, int W, double[][][] p,
			double[][][] d, double[][][] v, double[][][] r, double[] Q) throws IloException {
		// Model created
		final double V_GT = 0.98;
		final double V_ROT = 0.95;
		IloCplex cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1e-6);
		// cplex.setOut(null);

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
					objExpr = cplex.sum(objExpr, cplex.prod(x[s][k][t], r[s][k][t]));
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
//		IloNumExpr sum_xGT = cplex.constant(0);
//		IloNumExpr sum_dGT = cplex.constant(0);
//		IloNumExpr sum_xROT = cplex.constant(0);
//		IloNumExpr sum_dROT = cplex.constant(0);
//		for (int s = 0; s < S; s++) {
//			for (int k = 0; k < K; k++) {
//				for (int t = 0; t < T; t++) {
//					if (chunks.get(k).substring(0, 1).equals("R")) {
//						sum_xROT = cplex.sum(sum_xROT, x[s][k][t]);
//						sum_dROT = cplex.sum(sum_dROT, d[s][k][t]);
//
//					} else {
//						sum_xGT = cplex.sum(sum_xGT, x[s][k][t]);
//						sum_dGT = cplex.sum(sum_dGT, d[s][k][t]);
//
//					}
//
//				}
//
//			}
//
//		}
//		cplex.addGe(sum_xROT, cplex.prod(cplex.constant(V_ROT), sum_dROT));
//		cplex.addGe(sum_xGT, cplex.prod(cplex.constant(V_GT), sum_dGT));

		// Engaging solver and formatting console output
		cplex.solve();
		if (cplex.getStatus() == IloCplex.Status.Optimal)
			maxRelevance = cplex.getObjValue();
		else
			System.out.println("No optimal solution found in maxRelevance");
		cplex.close();
	}

}
