# siRNA_Scoring
This ReadMe file contains information on the data contained in the output files from the siRNA_Scoring.py script.


This file reports the average z-scores, MFEs, and EDs per nucleotide across a tiled user defined window of nucleotides. It reports which nucleotides are paired or unpaired across the tiled user defined window according to the ScanFold model, as well as the total number of base pairs in each tiled window. It reports the number of possible base pairs per nucleotide for all nucleotides in the tiled user defined window according the ScanFold-Scan results. It reports the reverse complement of the sequence in the tiled user defined window. It also reports the average of average z-scores, average MFEs, and average EDs for the tiled user defined window, as well as the range of nucleotides being analyzed in the tiled window 


The columns of data in the output file are as follows for each user defined window: 

1. pos1avgz
2. pos2avgz
3. pos3avgz	
4. pos4avgz 
5. pos5avgz	
6. pos6avgz	
7. pos7avgz	
8. pos8avgz
9. pos9avgz		
10. pos10avgz 	
11. pos11avgz 	
12. pos12avgz 	
13. pos13avgz 	
14. pos14avgz 	
15. pos15avgz 	
16. pos16avgz 	
17. pos17avgz 	
18. pos18avgz 	
19. pos19avgz 	
20. pos1avgMFE	
21. pos2avgMFE	
22. pos3avgMFE	
23. pos4avgMFE	
24. pos5avgMFE	
25. pos6avgMFE	
26. pos7avgMFE	
27. pos8avgMFE	
28. pos9avgMFE	
29. pos10avgMFE	
30. pos11avgMFE	
31. pos12avgMFE	
32. pos13avgMFE	
33. pos14avgMFE	
34. pos15avgMFE	
35. pos16avgMFE	
36. pos17avgMFE	
37. pos18avgMFE	
38. pos19avgMFE	
39. pos1avgED	
40. pos2avgED	
41. pos3avgED	
42. pos4avgED	
43. pos5avgED	
44. pos6avgED	
45. pos7avgED	
46. pos8avgED	
47. pos9avgED	
48. pos10avgED	
49. pos11avgED	
50. pos12avgED	
51. pos13avgED	
52. pos14avgED	
53. pos15avgED	
54. pos16avgED	
55. pos17avgED	
56. pos18avgED	
57. pos19avgED	
58. pos1pairedness	
59. pos2pairedness	
60. pos3pairedness	
61. pos4pairedness	
62. pos5pairedness	
63. pos6pairedness	
64. pos7pairedness	
65. pos8pairedness	
66. pos9pairedness	
67. pos10pairedness	
68. pos11pairedness	
69. pos12pairedness	
70. pos13pairedness	
71. pos14pairedness	
72. pos15pairedness	
73. pos16pairedness	
74. pos17pairedness	
75. pos18pairedness	
76. pos19pairedness	
77. #ofBPs	
78. pos1BP/nt	
79. pos2BP/nt	
80. pos3BP/nt	
81. pos4BP/nt	
82. pos5BP/nt	
83. pos6BP/nt	
84. pos7BP/nt	
85. pos8BP/nt	
86. pos9BP/nt	
87. pos10BP/nt	
88. pos11BP/nt	
89. pos12BP/nt	
90. pos13BP/nt	
91. pos14BP/nt	
92. pos15BP/nt	
93. pos16BP/nt	
94. pos17BP/nt	
95. pos18BP/nt	
96. pos19BP/nt	
97. Reverse Complement
98. zavg	
99. avgMFE	
100. avgED	
101. ntRange
102. TranscriptID


Below are the definition of each value type.

pos#avgzs: ScanFold average z-score value for each nucleotide position across the tiled user defined winodw size.

pos#avgMFE: Minumum free energy value for each nucleotide position across the tiled user defined winodw size.

pos#avgED: ScanFold average ensemble diversity value for each nucleotide position across the tiled user defined winodw size.

pos#pairedness: This tells whether each position across the tiled user defined window is paired or unpaired according to the ScanFold model. A value of 1 means that the nucleotide is paired and a value of 0 means the nucleotide is unpaired. 

#ofBP: The number of total number of base pairs in the tiled user defined window according the to ScanFold model.

pos#BP/nt: The number of possible base pairs each nucleotide in the user defined window was predicted to make according to ScanFold-Scan.

Reverse Complement: The reverse complement of the sequence in the tiled user defined window. Included for ease of identifying siRNA seed sequences.

zavg: Average of the average z-scores across the tile user defined window.

avgMFE: Average of the average minimum free energies across the tiled user defined window.

avgED: Average of the average ensemble diversities across the tiled user defined window.

ntRange: The range of nucleotides being analyzed in each window

TranscriptID: Transcript ID associated with each data window. This is included for ease of analysis when multiple files are concatenated into one larger file.
