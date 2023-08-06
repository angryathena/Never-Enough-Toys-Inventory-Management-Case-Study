import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;

import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;


public class assignment3 {

	private static double maxRevenue = 0;
	private static double maxRelevance = 0;
	private static double[][] yRev;
	private static double[][] yRel;
	private static double[][] xRev;
	private static double[][] xRel;
	private static double[][] gRev;
	private static double[][] gRel;
	private static double[][] IObjWeekly;
	private static double[][][] IObjWeeklyRaw;
	private static double[][] oldVolume;
	private static double[][][] oldVolumeRaw;
	private static int iter = 100;
	private static int mistakes = 0;
	private static int attempts = 0;

	private static ArrayList<ArrayList<Double>> simulate(ArrayList<String> chunks, double[][][] d18, double[][][] d19,
			double[][][] d, double[][][] p, double[][][] r, double[][][] v, int S, int K, int T, int W, double[][] y,
			double[] Q, double[][] gAgg) {
		ArrayList<ArrayList<Double>> solutions = new ArrayList<>();
		ArrayList<Double> simulatedRevenue = new ArrayList<>();
		ArrayList<Double> simulatedRelevance = new ArrayList<>();
		ArrayList<Double> simulatedGT = new ArrayList<>();
		ArrayList<Double> simulatedROT = new ArrayList<>();
		double[][][] dSim = new double[S][K][T];
		for (int i = 0; i < iter; i++) {
			// Random number generator and objective placeholders
			Random generator = new Random(i);
			double[][][] sales = new double[S][K][T];
			IObjWeeklyRaw = new double[S][K][T];
			// Simulating noisy demand
			for (int k = 0; k < K; k++) {
				for (int s = 0; s < S; s++) {
					// No uncertainty in the first week, orders = sales = demand, 0 leftover
					// inventory
					dSim[s][k][0] = d[s][k][0];
					for (int t = 1; t < T; t++) {
						// Computing standard deviation
						double avg = (d18[s][k][t] + d19[s][k][t]) / 2.0;
						double sigma = Math.sqrt((d18[s][k][t] - avg) * (d18[s][k][t] - avg)
								+ (d19[s][k][t] - avg) * (d19[s][k][t] - avg));
						// Adding noise
						double beta = generator.nextDouble();
						double gamma = generator.nextDouble();
						if (gamma < 0.5) {
							if (beta < 0.35)
								dSim[s][k][t] = d[s][k][t];
							if (beta >= 0.35 && beta < 0.65)
								dSim[s][k][t] = d[s][k][t] - sigma;
							if (beta >= 0.65 && beta < 0.95)
								dSim[s][k][t] = d[s][k][t] - 2 * sigma;
							if (beta >= 0.95)
								dSim[s][k][t] = d[s][k][t] - 3 * sigma;
						} else {
							if (beta < 0.35)
								dSim[s][k][t] = d[s][k][t];
							if (beta >= 0.35 && beta < 0.65)
								dSim[s][k][t] = d[s][k][t] + sigma;
							if (beta >= 0.65 && beta < 0.95)
								dSim[s][k][t] = d[s][k][t] + 2 * sigma;
							if (beta >= 0.95)
								dSim[s][k][t] = d[s][k][t] + 3 * sigma;
						}
						// Rounding demand for integrality
						dSim[s][k][t] = Math.max(0, Math.round(dSim[s][k][t]));
					}
				}
			}
			// Solving revenue weekly model
			double[][] oldSales = new double[S][K];
			double[][] oldInventory = new double[S][K];
			double[][][] g = new double[S][K][T];
			oldVolumeRaw = new double[S][K][W];
			for (int t = 0; t < T; t++) {
				for (int k = 0; k < K; k++) {
					for (int s = 0; s < S; s++) {
						if (t > 42)
							g[s][k][t] = gAgg[k][t-43];
						else
							g[s][k][t] = 1;
					}
				}
			}
			for (int t = 0; t < T; t++) {
				try {
					solveModel(chunks, t, S, K, T, W, p, d, v, g, oldSales, oldInventory, oldVolumeRaw, y, Q, 0);
					for (int s = 0; s < S; s++) {
						for (int k = 0; k < K; k++) {
							oldSales[s][k] = Math.min(IObjWeeklyRaw[s][k][t], dSim[s][k][t]);
							sales[s][k][t] = Math.min(IObjWeeklyRaw[s][k][t], dSim[s][k][t]);
							oldInventory[s][k] = IObjWeeklyRaw[s][k][t];
						}
					}

				} catch (IloException e) {
					System.out.println("A Cplex exception occured: " + e.getMessage());
					e.printStackTrace();
				}

			}
			// Computing performance metrics
			double totalRevenue = 0;
			double totalRelevance = 0;
			double salesGT = 0;
			double salesROT = 0;
			double demandGT = 0;
			double demandROT = 0;
			for (int s = 0; s < S; s++) {
				for (int k = 0; k < K; k++) {
					for (int t = 0; t < T; t++) {
						totalRevenue += sales[s][k][t] * p[s][k][t];
						totalRelevance += sales[s][k][t] * r[s][k][t];
						if (chunks.get(k).substring(0, 1).equals("R")) {
							salesROT += sales[s][k][t];
							demandROT += dSim[s][k][t];
						} else {
							salesGT += sales[s][k][t];
							demandGT += dSim[s][k][t];
						}
					}
				}
			}
			simulatedRevenue.add(totalRevenue);
			simulatedRelevance.add(totalRelevance);
			simulatedGT.add(salesGT / demandGT);
			simulatedROT.add(salesROT / demandROT);
		}
		solutions.add(simulatedRevenue);
		solutions.add(simulatedRelevance);
		solutions.add(simulatedGT);
		solutions.add(simulatedROT);
		return solutions;
	}

	public static void main(String[] args) throws FileNotFoundException {
		// Chunk list
		Scanner chunkList = new Scanner(new File("chunks2020.txt"));
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
		// Solving models for optimal sale schedules
		 revenueModel(chunks, 1, K, T, W, pAgg, dAgg, vAgg, rAgg, Q);
		 relevanceModel(chunks, 1, K, T, W, pAgg, dAgg, vAgg, rAgg, Q);

		// Reading previously-generated solutions to save time when simulating
//		yRev = new double[K][W];
//		gRev = new double[K][T];
//		yRel = new double[K][W];
//		gRel = new double[K][T];
//		Scanner gScan = new Scanner(new File("gRevenue.txt"));
//		for (int k = 0; k < K; k++) {
//			Scanner colGScan = new Scanner(gScan.nextLine());
//			for (int t = 0; t < 9; t++) {
//				gRev[k][t] = colGScan.nextDouble();
//			}
//			colGScan.close();
//		}
//		gScan = new Scanner(new File("gRelevance.txt"));
//		for (int k = 0; k < K; k++) {
//			Scanner colGScan = new Scanner(gScan.nextLine());
//			for (int t = 0; t < 9; t++) {
//				gRel[k][t] = colGScan.nextDouble();
//			}
//			colGScan.close();
//		}
//		gScan.close();
//		Scanner yScan = new Scanner(new File("yRevenue.txt"));
//		for (int k = 0; k < K; k++) {
//			Scanner colYScan = new Scanner(yScan.nextLine());
//			for (int w = 0; w < W; w++) {
//				yRev[k][w] = colYScan.nextDouble();
//			}
//			colYScan.close();
//		}
//		yScan = new Scanner(new File("yRelevance.txt"));
//		for (int k = 0; k < K; k++) {
//			Scanner colYScan = new Scanner(yScan.nextLine());
//			for (int w = 0; w < W; w++) {
//				yRel[k][w] = colYScan.nextDouble();
//			}
//			colYScan.close();
//		}
//		yScan.close();

		// Print relevant large-model output for simulation use
			PrintStream ps1 = new PrintStream(new File("yRevenue.txt"));
			PrintStream ps2 = new PrintStream(new File("yRelevance.txt"));
			PrintStream ps3 = new PrintStream(new File("gRevenue.txt"));
			PrintStream ps4 = new PrintStream(new File("gRelevance.txt"));
				for (int k=0; k<K; k++) {
					for (int t=0; t<9; t++) {
						ps3.print(gRev[k][t] + " ");
						ps4.print(gRel[k][t] + " ");
					}
					ps3.println();
					ps4.println();
				}
			for (int k=0; k<K; k++) {
				for (int w=0; w<W; w++) {
					ps1.print(yRev[k][w] + " ");
					ps2.print(yRel[k][w] + " ");
				}
				ps1.println();
				ps2.println();
			}
			ps1.close();
			ps2.close();
			ps3.close();
			ps4.close();

		// Simulating for REVENUE-maximizing optimal solution
		ArrayList<ArrayList<Double>> simulated_xRev = simulate(chunks, d18, d19, d20, p20, r20, v20, S, K, T, W, yRev,
				Q, gRev);
		// Simulating for RELEVANCE-maximizing optimal solution
		ArrayList<ArrayList<Double>> simulated_xRel = simulate(chunks, d18, d19, d20, p20, r20, v20, S, K, T, W, yRel,
				Q, gRel);

		// Printing results
		PrintStream ps5 = new PrintStream(new File("SIM_xRev.txt"));
		ps5.println("Revenue" + "	" + "Relevance" + "	" + "Service GT" + "	" + "Service ROT");
		for (int i = 0; i < iter; i++)
			ps5.println(simulated_xRev.get(0).get(i) + "	" + simulated_xRev.get(1).get(i) + "	"
					+ simulated_xRev.get(2).get(i) + "	" + simulated_xRev.get(3).get(i));
		ps5.close();
		PrintStream ps6 = new PrintStream(new File("SIM_xRel.txt"));
		ps6.println("Revenue" + "	" + "Relevance" + "	" + "Service GT" + "	" + "Service ROT");
		for (int i = 0; i < iter; i++)
			ps6.println(simulated_xRel.get(0).get(i) + " " + simulated_xRel.get(1).get(i) + " "
					+ simulated_xRel.get(2).get(i) + " " + simulated_xRel.get(3).get(i));
		ps6.close();

		// Console: total revenue and relevance of noiseless model
//			System.out.println("--------------------------------------------------------------------------");
//			System.out.println("Total Revenue: " + maxRevenue);
//			System.out.println("Total Relevance: " + maxRelevance);
		} catch (IloException e) {
			System.out.println("A Cplex exception occured: " + e.getMessage());
			e.printStackTrace();
		}

		System.out.println(mistakes);
		System.out.println(attempts);

	}

	public static void revenueModel(ArrayList<String> chunks, int S, int K, int T, int W, double[][] pAgg,
			double[][] dAgg, double[][] vAgg, double[][] rAgg, double[] Q) throws IloException {
		// Model created
		final double V_GT = 0.98;
		final double V_ROT = 0.95;
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
					sum_xv[k][t] = cplex.sum(sum_xv[k][t], cplex.prod(cplex.sum(I[k][t - 43], z[k][t-43]), vAgg[k][t]));
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

		// Service Level Constraints: 4j
//		IloNumExpr sum_xGT = cplex.constant(0);
//		IloNumExpr sum_dGT = cplex.constant(0);
//		IloNumExpr sum_xROT = cplex.constant(0);
//		IloNumExpr sum_dROT = cplex.constant(0);
//		for (int s = 0; s < S; s++) {
//			for (int k = 0; k < K; k++) {
//				for (int t = 0; t < T; t++) {
//					if (chunks.get(k).substring(0, 1).equals("R")) {
//						if (t>=43)
//							sum_xROT = cplex.sum(sum_xROT, z[k][t-43]);
//						else
//							sum_xROT = cplex.sum(sum_xROT, x[k][t]);
//						sum_dROT = cplex.sum(sum_dROT, dAgg[k][t]);
//
//					} else {
//						if (t>=43)
//						sum_xGT = cplex.sum(sum_xGT, z[k][t-43]);
//						else
//							sum_xGT = cplex.sum(sum_xGT, x[k][t]);
//						sum_dGT = cplex.sum(sum_dGT, dAgg[k][t]);
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
		if (cplex.getStatus() == IloCplex.Status.Optimal) {
			maxRevenue = cplex.getObjValue();
			yRev = new double[K][W];
			gRev = new double[K][9];
			xRev = new double[K][T];
			for (int k = 0; k < K; k++) {
				for (int t = 0; t < 9; t++) {
					gRev[k][t] = cplex.getValue(g[k][t]);
				}
			}
			for (int k = 0; k < K; k++) {
				for (int t = 0; t < T; t++) {
					xRev[k][t] = cplex.getValue(x[k][t]);
				}
			}
			for (int w = 0; w < W; w++) {
				for (int k = 0; k < K; k++) {
					yRev[k][w] = cplex.getValue(y[k][w]);
				}
			}
		} else
			System.out.println("No optimal solution found in maxRevenue");
		cplex.close();
	}

	public static void relevanceModel(ArrayList<String> chunks, int S, int K, int T, int W, double[][] pAgg,
			double[][] dAgg, double[][] vAgg, double[][] rAgg, double[] Q) throws IloException {
		// Model created
		final double V_GT = 0.98;
		final double V_ROT = 0.95;
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
					objExpr = cplex.sum(objExpr, cplex.prod(x[k][t], rAgg[k][t]));
					cplex.addLe(x[k][t], dAgg[k][t]);
				}
				for (int t = 43; t < T; t++) {
					objExpr = cplex.sum(objExpr, cplex.prod(z[k][t - 43], rAgg[k][t]));
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
					sum_xv[k][t] = cplex.sum(sum_xv[k][t], cplex.prod(cplex.sum(I[k][t - 43], z[k][t-43]), vAgg[k][t]));
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

		// Service Level Constraints: 4j
//		IloNumExpr sum_xGT = cplex.constant(0);
//		IloNumExpr sum_dGT = cplex.constant(0);
//		IloNumExpr sum_xROT = cplex.constant(0);
//		IloNumExpr sum_dROT = cplex.constant(0);
//		for (int s = 0; s < S; s++) {
//			for (int k = 0; k < K; k++) {
//				for (int t = 0; t < T; t++) {
//					if (chunks.get(k).substring(0, 1).equals("R")) {
//						if (t>=43)
//							sum_xROT = cplex.sum(sum_xROT, z[k][t-43]);
//						else
//							sum_xROT = cplex.sum(sum_xROT, x[k][t]);
//						sum_dROT = cplex.sum(sum_dROT, dAgg[k][t]);
//
//					} else {
//						if (t>=43)
//						sum_xGT = cplex.sum(sum_xGT, z[k][t-43]);
//						else
//							sum_xGT = cplex.sum(sum_xGT, x[k][t]);
//						sum_dGT = cplex.sum(sum_dGT, dAgg[k][t]);
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
		if (cplex.getStatus() == IloCplex.Status.Optimal) {
			maxRelevance = cplex.getObjValue();
			yRel = new double[K][W];
			gRel = new double[K][9];
			xRel = new double[K][T];
			for (int k = 0; k < K; k++) {
				for (int t = 0; t < 9; t++) {
					gRel[k][t] = cplex.getValue(g[k][t]);
				}
			}
			for (int k = 0; k < K; k++) {
				for (int t = 0; t < T; t++) {
					xRel[k][t] = cplex.getValue(x[k][t]);
				}
			}
			for (int w = 0; w < W; w++) {
				for (int k = 0; k < K; k++) {
					yRel[k][w] = cplex.getValue(y[k][w]);
				}
			}
		} else
			System.out.println("No optimal solution found in maxRevenue");
		cplex.close();
	}

	public static void solveModel(ArrayList<String> chunks, int t, int S, int K, int T, int W, double[][][] param,
			double[][][] d, double[][][] v, double[][][] g, double[][] zLast, double[][] ILast, double qLast[][][],
			double[][] y, double[] Q, int safe) throws IloException {
		// Model created
		IloCplex cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1e-3);
		// cplex.setOut(null);

		// Decision variables
		IloNumVar[][] x = new IloNumVar[S][K]; // size
		for (int k = 0; k < K; k++) {
			for (int s = 0; s < S; s++) {
				x[s][k] = cplex.intVar(0, Integer.MAX_VALUE);
			}
		}

		// Objective function
		IloNumExpr[][] I = new IloNumExpr[S][K];

		for (int k = 0; k < K; k++) {
			for (int s = 0; s < S; s++) {
				I[s][k] = cplex.sum(ILast[s][k] - zLast[s][k], x[s][k]);
			}
		}

		IloNumExpr objExpr = cplex.constant(0);
		for (int k = 0; k < K; k++) {
			for (int s = 0; s < S; s++) {
				objExpr = cplex.sum(objExpr, cplex.prod(I[s][k], param[s][k][t]));
			}
		}
		cplex.addMaximize(objExpr);

		// Upper bounds on order size (no safety stock yet)

		for (int k = 0; k < K; k++) {
			for (int s = 0; s < S; s++) {
				if (t >= 43 && t <= 50)
					cplex.addLe(x[s][k],
							Math.max(0, g[s][k][t - 43] * (d[s][k][t] + d[s][k][t + 1] - ILast[s][k] + zLast[s][k])));
				else if (t == 51)
					cplex.addLe(x[s][k], Math.max(0, g[s][k][t - 43] * (d[s][k][t] - ILast[s][k] + zLast[s][k])));
				else
					cplex.addLe(x[s][k], Math.max(0, d[s][k][t] - ILast[s][k] + zLast[s][k]));
			}
		}

		// Capacity constraints
		IloNumExpr[][][] q = new IloNumExpr[S][K][W];
		IloNumExpr[] qSum = new IloNumExpr[W];
		for (int w = 0; w < W; w++) {
			qSum[w] = cplex.constant(0);
			for (int k = 0; k < K; k++) {
				for (int s = 0; s < S; s++) {
					if (t > 0)
						q[s][k][w] = cplex.sum(qLast[s][k][w], cplex.diff(cplex.prod(x[s][k], v[s][k][t] * y[k][w]),
								zLast[s][k] * v[s][k][t - 1] * y[k][w]));
					else
						q[s][k][w] = cplex.prod(x[s][k], v[s][k][t] * y[k][w]);
					qSum[w] = cplex.sum(qSum[w], q[s][k][w]);
				}
			}
			cplex.addLe(qSum[w], Q[w]);
		}

		// Engaging solver and formatting console output
		cplex.solve();
		attempts++;
		// cplex.exportModel("Model.lp");
		if (cplex.getStatus() == IloCplex.Status.Optimal) {
			for (int k = 0; k < K; k++) {
				for (int s = 0; s < S; s++) {
					IObjWeeklyRaw[s][k][t] = cplex.getValue(I[s][k]);
					for (int w = 0; w < W; w++) {
						oldVolumeRaw[s][k][w] = cplex.getValue(q[s][k][w]);
					}
				}
			}
		} else {
			System.out.println("No optimal solution found in maxRevenue");
			mistakes++;
		}
		cplex.close();
	}
}
