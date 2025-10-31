================================================================================
                    PARAMETRIC STUDY FIGURES - GUIDE
                          Heat Exchanger Optimization
================================================================================

Generated: 28-Oct-2025
Total Studies: 7
Figures per Study: 6 subplots
File Formats: PNG (images), FIG (editable MATLAB), MAT (raw data)

TEMPERATURES CORRECTED:
  Hot fluid inlet:  20.00°C (293.15 K) - EXACTLY
  Cold fluid inlet:  8.00°C (281.15 K) - EXACTLY

================================================================================
                           FIGURE DESCRIPTIONS
================================================================================

Each parametric study figure contains 6 subplots showing the effect of
varying one design parameter:

1. Overall Heat Transfer Coefficient (U) [W/m²K]
   → Higher is better (more efficient heat transfer)
   → Target: 80-100 W/m²K

2. Total Heat Transfer Area [m²]
   → Lower is better (more compact, less material cost)
   → Target: Minimize while meeting heat duty

3. Total Number of Tubes
   → Affects manufacturability and cost
   → Result of area requirement ÷ area per tube

4. Tube-Side Pressure Drop [kPa]
   → Lower is better (less pumping power for water)
   → Target: < 5-10 kPa (acceptable for water systems)

5. Shell-Side Pressure Drop [kPa]
   → Lower is better (less fan power for air)
   → Target: < 100 kPa (acceptable for air systems)

6. Shell Diameter [m]
   → Affects overall size and cost
   → Smaller is generally better for compactness

================================================================================
                        STUDY 1: BAFFLE SPACING RATIO
================================================================================

Files: 1_bafflespacingratio.png, .fig, .mat
Variable Range: 0.3 to 0.8 (dimensionless)
Current Optimal: 0.8

DESCRIPTION:
Baffle spacing as a fraction of shell diameter. Baffles direct the shell-side
flow across the tubes to improve heat transfer.

KEY FINDINGS:
- U remains relatively constant (~89-90 W/m²K) - baffle spacing has minimal
  effect on overall heat transfer coefficient
- Shell-side pressure drop DECREASES DRAMATICALLY as baffle spacing increases
  • At 0.3: ~1100 kPa (UNACCEPTABLE - would require huge fan power)
  • At 0.8: ~53 kPa (ACCEPTABLE - reasonable fan power)
- Tube-side pressure drop remains constant (baffles don't affect tube side)

RECOMMENDATION:
Use 0.8 (wider baffle spacing) to minimize shell-side pressure drop while
maintaining good heat transfer. This was the CRITICAL fix from the original
design which had 0.4 ratio causing 2177 kPa pressure drop!

================================================================================
                          STUDY 2: FIN PITCH
================================================================================

Files: 2_pf.png, .fig, .mat
Variable Range: 0.003 to 0.006 m (3mm to 6mm)
Current Optimal: 0.005 m (5mm)

DESCRIPTION:
Distance between fins on the tubes. Smaller pitch = more fins per meter.

FIN DENSITY:
- 3mm pitch: 333 fins/meter (very dense)
- 5mm pitch: 200 fins/meter (optimal)
- 6mm pitch: 167 fins/meter (sparse)

KEY FINDINGS:
- U INCREASES as fin pitch increases (counterintuitive!)
  • At 3mm: ~72 W/m²K
  • At 5mm: ~86 W/m²K
  • At 6mm: ~90 W/m²K
- Total area DECREASES (fewer but more effective fins)
- Shell-side pressure drop DECREASES significantly
  • At 3mm: ~84 kPa
  • At 5mm: ~57 kPa
  • At 6mm: ~50 kPa

EXPLANATION:
Too many fins causes flow resistance and reduced fin efficiency. Optimal
spacing allows better air flow between fins.

RECOMMENDATION:
5mm pitch provides best balance - good U coefficient with acceptable pressure
drop. Going to 6mm gives marginal improvements.

================================================================================
                     STUDY 3: NUMBER OF TUBE PASSES
================================================================================

Files: 3_Np.png, .fig, .mat
Variable Range: 2, 4, 6, 8 passes
Current Optimal: 4 passes

DESCRIPTION:
Number of times the tube-side fluid travels the length of the exchanger. More
passes = higher velocity in each tube.

KEY FINDINGS:
- U INCREASES strongly with more passes
  • 2 passes: ~61 W/m²K (poor)
  • 4 passes: ~86 W/m²K (good)
  • 8 passes: ~103 W/m²K (excellent)
- Total area decreases (more efficient heat transfer)
- Tube-side pressure drop INCREASES DRAMATICALLY
  • 2 passes: 0.08 kPa (very low)
  • 4 passes: 0.99 kPa (low)
  • 8 passes: 10.1 kPa (moderate)
- Shell-side pressure drop also increases (smaller bundle, tighter spacing)

RECOMMENDATION:
4 passes provides excellent balance. Going to 6 or 8 passes gives marginal
U improvement but quadruples pressure drop.

================================================================================
                     STUDY 4: TUBE INNER DIAMETER
================================================================================

Files: 4_Dinternal.png, .fig, .mat
Variable Range: 0.019 to 0.022 m (19mm to 22mm ID)
Current Optimal: 0.022 m (22mm) - STANDARD SIZE

DESCRIPTION:
Inner diameter of tubes. All values tested are STANDARD sizes for 25mm OD
tube with different wall thicknesses.

STANDARD SIZES TESTED:
- 19mm ID: 25mm OD x 3mm wall (heavy wall)
- 20mm ID: 25mm OD x 2.5mm wall
- 21mm ID: 25mm OD x 2mm wall
- 22mm ID: 25mm OD x 1.5mm wall (thin wall)

KEY FINDINGS:
- U DECREASES slightly as diameter increases
  • 19mm: ~86 W/m²K
  • 22mm: ~83 W/m²K
- Tube-side pressure drop DECREASES significantly
  • 19mm: ~0.99 kPa
  • 22mm: ~0.49 kPa
- All pressure drops are very low and acceptable

RECOMMENDATION:
22mm ID (1.5mm wall) provides lowest pressure drop with minimal sacrifice
in U. All tested sizes perform well - choose based on availability and cost.

================================================================================
                         STUDY 5: TUBE LENGTH
================================================================================

Files: 5_Ltube.png, .fig, .mat
Variable Range: 1.0 to 2.0 m
Current Optimal: 1.0 m

DESCRIPTION:
Length of individual tubes. Longer tubes = more area per tube = fewer tubes
needed, but also longer overall exchanger.

KEY FINDINGS:
- U DECREASES as length increases
  • 1.0m: ~91 W/m²K (best)
  • 2.0m: ~82 W/m²K
- Number of tubes DECREASES (longer tubes provide more area each)
  • 1.0m: 61 tubes
  • 2.0m: 34 tubes
- Tube-side pressure drop INCREASES (longer flow path)
  • 1.0m: ~1.0 kPa
  • 2.0m: ~4.6 kPa
- Shell-side pressure drop INCREASES DRAMATICALLY
  • 1.0m: ~78 kPa
  • 2.0m: ~430 kPa

RECOMMENDATION:
1.0m tubes provide best performance and acceptable pressure drops. Shorter
exchanger is also easier to fit in mechanical rooms.

================================================================================
                          STUDY 6: FIN HEIGHT
================================================================================

Files: 6_lf.png, .fig, .mat
Variable Range: 0.003 to 0.008 m (3mm to 8mm)
Current Optimal: 0.005 m (5mm)

DESCRIPTION:
Height of fins extending from tube outer diameter. Taller fins = more surface
area but potentially lower fin efficiency.

KEY FINDINGS:
- U DECREASES slightly as fin height increases
  • 3mm: ~94 W/m²K
  • 5mm: ~89 W/m²K
  • 8mm: ~83 W/m²K
- Total area INCREASES (more fin surface area)
- Shell-side pressure drop varies moderately
- Diminishing returns due to reduced fin efficiency for tall fins

RECOMMENDATION:
5mm height provides good balance. Taller fins add area but become less
efficient at transferring heat from their tips.

================================================================================
                        STUDY 7: FIN THICKNESS
================================================================================

Files: 7_tf.png, .fig, .mat
Variable Range: 0.0005 to 0.002 m (0.5mm to 2.0mm)
Current Optimal: 0.001 m (1.0mm)

DESCRIPTION:
Thickness of fin material. Thicker fins = better heat conduction along fin
but fewer fins can fit in same space.

KEY FINDINGS:
- U INCREASES with fin thickness
  • 0.5mm: ~87 W/m²K
  • 1.0mm: ~89 W/m²K
  • 2.0mm: ~92 W/m²K
- Changes are relatively moderate across the range
- All thicknesses show good performance

RECOMMENDATION:
1.0mm provides good fin efficiency and is standard for aluminum fins.
Thicker fins give marginal improvement but use more material.

================================================================================
                         HOW TO USE THESE FIGURES
================================================================================

FOR REPORT WRITING:
1. Use PNG files in your Word/LaTeX documents
2. Refer to figures to justify design choices
3. Explain trade-offs (e.g., "4 passes chosen to balance U vs pressure drop")

FOR PRESENTATIONS:
1. PNG files work in PowerPoint
2. Focus on 1-2 key subplots per study
3. Use to show optimization process

FOR FURTHER ANALYSIS:
1. Load MAT files in MATLAB: load('1_bafflespacingratio.mat')
2. Edit FIG files to customize appearance
3. Extract specific data for tables

FOR DESIGN REFINEMENT:
If constraints change (e.g., pressure drop limits), review figures to find
alternative optimal points.

================================================================================
                            QUICK REFERENCE
================================================================================

OPTIMIZATION SEQUENCE USED:
1. Fixed baffle spacing (0.4 → 0.8) - reduced ΔP_s from 2177 → 50 kPa
2. Optimized fin pitch (3mm → 5mm) - increased U from 72 → 86 W/m²K
3. Optimized tube passes (6 → 4) - balanced U vs ΔP_t
4. Selected standard tube size (22mm ID) - lowest ΔP_t
5. Optimized tube length (1.5m → 1.0m) - best U, lower ΔP_s
6. Verified fin dimensions (5mm height, 1mm thickness)
7. Applied LMTD correction factor (F = 0.96)

FINAL DESIGN PARAMETERS:
- Tube: 25mm OD x 22mm ID x 1.0m (EU Standard)
- Fins: 5mm height, 1mm thick, 5mm pitch
- Passes: 4
- Baffle spacing ratio: 0.8
- LMTD correction: F = 0.96

FINAL PERFORMANCE:
- U = 88.1 W/m²K
- Area = 17.6 m²
- ΔP_tube = 0.47 kPa ✓
- ΔP_shell = 50 kPa ✓
- 65 tubes total

================================================================================
                              END OF GUIDE
================================================================================
