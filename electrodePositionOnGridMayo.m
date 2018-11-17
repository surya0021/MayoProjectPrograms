function [row,column,electrodeArray] = electrodePositionOnGridMayo(electrodeNum,monkeyName)

if strcmp(monkeyName,'w')
 electrodeArray = ...
   [46    50    15    17     9     6   0  75    73    71    69    67    65
    11     8    90    89    55    56   0  72    70    68    66    79    77
    57    58    52    54    19    13   0   2    81    80    78    76    74
    25    23    12    10    92    91   0  41    39    37    35    33     1
    93    94    59    60    62    21   0  34     3    83    84    82    45
    31    29    27    20    16    14   0  49    47    43    40    38    36
    22    18    96    63    61    64   0  48     5     7     4    85    86
    95    32    30    28    26    24   0  88    87    53    51    44    42];
elseif strcmp(monkeyName,'a')
  electrodeArray = ...
   [75    73    71    69    67    65   0  46    50    15    17     9     6 
    72    70    68    66    79    77   0  11     8    90    89    55    56
     2    81    80    78    76    74   0  57    58    52    54    19    13
    41    39    37    35    33     1   0  25    23    12    10    92    91 
    34     3    83    84    82    45   0  93    94    59    60    62    21
    49    47    43    40    38    36   0  31    29    27    20    16    14 
    48     5     7     4    85    86   0  22    18    96    63    61    64 
    88    87    53    51    44    42   0  95    32    30    28    26    24]; 
end
if electrodeNum<1 || electrodeNum>96
    disp('Electrode Number out of range');
    row=1;column=1;
else
    [row,column] = find(electrodeArray==electrodeNum);
end
end