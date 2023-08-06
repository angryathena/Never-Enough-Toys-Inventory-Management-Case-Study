import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Scanner;

import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.MIPStartEffort;

public class assignment1 {

	public static double revenueEq = 0;
	public static double revenueLe = 0;
	public static double totalSalesROT = 0;
	public static double totalDemandROT = 0;
	public static double totalSalesGT = 0;
	public static double totalDemandGT = 0;
	public static int[][] yLe = new int[495][2];
	public static int[][][] xLe = new int[9][495][52];
	public static double[][][] qLe = new double[495][2][52];
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
				}
				colVolume.close();
			}
		}
		volume.close();

		// Chunk list 2019
		Scanner chunkList2019 = new Scanner(new File("chunks2019.txt"));
		ArrayList<String> chunks2019 = new ArrayList<>();
		while (chunkList2019.hasNextLine()) {
			chunks2019.add(chunkList2019.nextLine());
		}
		chunkList2019.close();

		// Demand list 2019
		Scanner demand2019 = new Scanner(new File("demands2019a1.txt"));
		double[][][] d2019 = new double[S][K][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				Scanner colDemand = new Scanner(demand2019.nextLine());
				for (int t = 0; t < T; t++) {
					d2019[s][k][t] = colDemand.nextDouble();
				}
				colDemand.close();
			}
		}
		demand2019.close();

		// Revenue list 2019
		Scanner profit2019 = new Scanner(new File("revenues2019a1.txt"));
		double[][][] p2019 = new double[S][K][T];
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				Scanner colProfit = new Scanner(profit2019.nextLine());
				for (int t = 0; t < T; t++) {
					p2019[s][k][t] = colProfit.nextDouble();
				}
				colProfit.close();
			}
		}
		profit2019.close();

		double[] Q = { 2250, 450 };

		try {

			solveModel(chunks, S, K, T, W, p, d, v, Q);
			// testingModel(xEq, yEq, qEq, S, K, T, W, p, d, v, Q);

		} catch (IloException e) {
			System.out.println("A Cplex exception occured: " + e.getMessage());
			e.printStackTrace();
		}

		try {
			PrintStream ps1 = new PrintStream(new File("allSales.txt"));
			for (int s = 0; s < S; s++) {
				for (int k = 0; k < K; k++) {
					for (int t = 0; t < T; t++) {
						ps1.print(xEq[s][k][t] + " ");
					}
					ps1.println();
				}
			}
			ps1.close();

			PrintStream ps2 = new PrintStream(new File("warehouseChoice.txt"));
			for (int k = 0; k < K; k++) {
				if (yEq[k][0] == 1)
					ps2.println(0);
				else
					ps2.println(1);
			}
			ps2.close();

			PrintStream ps3 = new PrintStream(new File("occupationLevel.txt"));
			for (int t = 0; t < T; t++) {
				double totalVolumeA = 0;
				double totalVolumeB = 0;
				for (int k = 0; k < K; k++) {
					totalVolumeA += qEq[k][0][t];
					totalVolumeB += qEq[k][1][t];
				}
				ps3.println(totalVolumeA / Q[0] + " " + totalVolumeB / Q[1]);
			}
			ps3.close();

			PrintStream ps4 = new PrintStream(new File("serviceLevel.txt"));
			for (int t = 0; t < T; t++) {
				double weeklySalesROT = 0;
				double weeklySalesGT = 0;
				double weeklyDemandROT = 0;
				double weeklyDemandGT = 0;
				for (int k = 0; k < K; k++) {
					for (int s = 0; s < S; s++) {
						if (chunks.get(k).substring(0, 1).equals("R")) {
							weeklySalesROT += xEq[s][k][t];
							weeklyDemandROT += d[s][k][t];

						} else {
							weeklySalesGT += xEq[s][k][t];
							weeklyDemandGT += d[s][k][t];

						}

					}

				}
				ps4.println((double) weeklySalesGT / weeklyDemandGT + " " + (double) weeklySalesROT / weeklyDemandROT);
			}

			// 2019 stuff
			double rev19 = 0;
			double salesGT19 = 0;
			double salesROT19 = 0;
			double demandGT19 = 0;
			double demandROT19 = 0;
			for (int k = 0; k < K; k++) {
				for (int s = 0; s < S; s++) {
					for (int t = 0; t < T; t++) {
						rev19 += Math.min(xEq[s][k][t], d2019[s][k][t]) * p2019[s][k][t];
						if (chunks.get(k).substring(0, 1).equals("R")) {
							salesROT19 += Math.min(xEq[s][k][t], d2019[s][k][t]);
							demandROT19 += d2019[s][k][t];
						} else {
							salesGT19 += Math.min(xEq[s][k][t], d2019[s][k][t]);
							demandGT19 += d2019[s][k][t];
						}
					}
				}

			}
			System.out.println("2019 Revenue: " + rev19);
			System.out.println("2019 Service Level: " + 100 * salesGT19 / demandGT19 + "% GT, "
					+ 100 * salesROT19 / demandROT19 + "% ROT");

		} catch (IOException e) {
			System.out.println("An error occurred on printing.");
			e.printStackTrace();
		}

	}

	public static void solveModel(ArrayList<String> chunks, int S, int K, int T, int W, double[][][] p, double[][][] d,
			double[][][] v, double[] Q) throws IloException {
		// Model created
		IloCplex cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1e-6);
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

		// Engaging solver and formatting console output
		cplex.solve();
		if (cplex.getStatus() == IloCplex.Status.Optimal) {
			System.out.println("Found optimal solution!");
			System.out.println("Objective = " + cplex.getObjValue());
			revenueEq = cplex.getObjValue();

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
			System.out.println("Service Level: " + 100 * totalSalesGT / totalDemandGT + "% GT, "
					+ 100 * totalSalesROT / totalDemandROT + "% ROT");
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

	public static void testingModel(int[][][] xVals, int[][] yVals, double[][][] qVals, int S, int K, int T, int W,
			double[][][] pi, double[][][] d, double[][][] v, double[] Q) throws IloException {
		// Model created
		IloCplex cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1e-6);

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
					objExpr = cplex.sum(objExpr, cplex.prod(x[s][k][t], pi[s][k][t]));
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
				cplex.addLe(sum_xv[k][t], sum_w_q[k][t]);
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

		// warm start

		IloNumVar[] vars = new IloNumVar[S * K * T + K * W + K * W * T];
		double[] vals = new double[S * K * T + K * W + K * W * T];
		int index = 0;
		for (int s = 0; s < S; s++) {
			for (int k = 0; k < K; k++) {
				for (int t = 0; t < T; t++) {
					vars[index] = x[s][k][t];
					vals[index] = xVals[s][k][t];
					index++;
				}

			}
		}
		for (int k = 0; k < K; k++) {
			for (int w = 0; w < W; w++) {
				vars[index] = y[k][w];
				vals[index] = yVals[k][w];
				index++;
			}

		}
		for (int k = 0; k < K; k++) {
			for (int w = 0; w < W; w++) {
				for (int t = 0; t < T; t++) {
					vars[index] = q[k][w][t];
					vals[index] = qVals[k][w][t];
					index++;
				}
			}
		}

		cplex.setStart(vals, vals, vars, vals, vals, null);
		MIPStartEffort effort = IloCplex.MIPStartEffort.SolveFixed;
		cplex.addMIPStart(vars, vals, effort);

		// Engaging solver and formatting console output

		cplex.solve();
		if (cplex.getStatus() == IloCplex.Status.Optimal) {
			System.out.println("Found optimal solution!");
			System.out.println("Objective = " + cplex.getObjValue());
			revenueLe = cplex.getObjValue();
			System.out.println("Equality Objective: " + revenueEq);
			double totalSales = 0;
			double totalDemand = 0;
			for (int s = 0; s < S; s++) {
				for (int k = 0; k < K; k++) {
					for (int t = 0; t < T; t++) {
						totalSales += cplex.getValue(x[s][k][t]);
						totalDemand += d[s][k][t];
						xLe[s][k][t] = (int) cplex.getValue(x[s][k][t]);
					}
				}
			}
			System.out.println("Service Level: " + totalSales / totalDemand + "%");
			for (int k = 0; k < K; k++) {
				for (int w = 0; w < W; w++) {
					yLe[k][w] = (int) cplex.getValue(y[k][w]);
				}
			}
			for (int k = 0; k < K; k++) {
				for (int w = 0; w < W; w++) {
					for (int t = 0; t < T; t++) {
						qLe[k][w][t] = cplex.getValue(q[k][w][t]);
					}

				}
			}

		} else {
			System.out.println("No optimal solution found");
		}
		cplex.close();
	}

}
