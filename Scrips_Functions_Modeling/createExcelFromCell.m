function createExcelFromCell(A_cell, filename)

alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

% Write the arrays to the Excel file, each array in a separate column
for i = 1:length(A_cell)

    if i > 26*26
        range = [ 'Z' alphabet(i-26*26) num2str(1) ':Z' alphabet(i-26*26) num2str(length(A_cell{i})) ];
    elseif i > 26*25
        range = [ 'Y' alphabet(i-26*25) num2str(1) ':Y' alphabet(i-26*25) num2str(length(A_cell{i})) ];
    elseif i > 26*24
        range = [ 'X' alphabet(i-26*24) num2str(1) ':X' alphabet(i-26*24) num2str(length(A_cell{i})) ];
    elseif i > 26*23
        range = [ 'W' alphabet(i-26*23) num2str(1) ':W' alphabet(i-26*23) num2str(length(A_cell{i})) ];
    elseif i > 26*22
        range = [ 'V' alphabet(i-26*22) num2str(1) ':V' alphabet(i-26*22) num2str(length(A_cell{i})) ];
    elseif i > 26*21
        range = [ 'U' alphabet(i-26*21) num2str(1) ':U' alphabet(i-26*21) num2str(length(A_cell{i})) ];
    elseif i > 26*20
        range = [ 'T' alphabet(i-26*20) num2str(1) ':T' alphabet(i-26*20) num2str(length(A_cell{i})) ];
    elseif i > 26*19
        range = [ 'S' alphabet(i-26*19) num2str(1) ':S' alphabet(i-26*19) num2str(length(A_cell{i})) ];
    elseif i > 26*18
        range = [ 'R' alphabet(i-26*18) num2str(1) ':R' alphabet(i-26*18) num2str(length(A_cell{i})) ];
    elseif i > 26*17
        range = [ 'Q' alphabet(i-26*17) num2str(1) ':Q' alphabet(i-26*17) num2str(length(A_cell{i})) ];
    elseif i > 26*16
        range = [ 'P' alphabet(i-26*16) num2str(1) ':P' alphabet(i-26*16) num2str(length(A_cell{i})) ];
    elseif i > 26*15
        range = [ 'O' alphabet(i-26*15) num2str(1) ':O' alphabet(i-26*15) num2str(length(A_cell{i})) ];
    elseif i > 26*14
        range = [ 'N' alphabet(i-26*14) num2str(1) ':N' alphabet(i-26*14) num2str(length(A_cell{i})) ];
    elseif i > 26*13
        range = [ 'M' alphabet(i-26*13) num2str(1) ':M' alphabet(i-26*13) num2str(length(A_cell{i})) ];
    elseif i > 26*12
        range = [ 'L' alphabet(i-26*12) num2str(1) ':L' alphabet(i-26*12) num2str(length(A_cell{i})) ];
    elseif i > 26*11
        range = [ 'K' alphabet(i-26*11) num2str(1) ':K' alphabet(i-26*11) num2str(length(A_cell{i})) ];
    elseif i > 26*10
        range = [ 'J' alphabet(i-26*10) num2str(1) ':J' alphabet(i-26*10) num2str(length(A_cell{i})) ];
    elseif i > 26*9
        range = [ 'I' alphabet(i-26*9) num2str(1) ':I' alphabet(i-26*9) num2str(length(A_cell{i})) ];
    elseif i > 26*8
        range = [ 'H' alphabet(i-26*8) num2str(1) ':H' alphabet(i-26*8) num2str(length(A_cell{i})) ];
    elseif i > 26*7
        range = [ 'G' alphabet(i-26*7) num2str(1) ':G' alphabet(i-26*7) num2str(length(A_cell{i})) ];
    elseif i > 26*6
        range = [ 'F' alphabet(i-26*6) num2str(1) ':F' alphabet(i-26*6) num2str(length(A_cell{i})) ];
    elseif i > 26*5
        range = [ 'E' alphabet(i-26*5) num2str(1) ':E' alphabet(i-26*5) num2str(length(A_cell{i})) ];
    elseif i > 26*4
        range = [ 'D' alphabet(i-26*4) num2str(1) ':D' alphabet(i-26*4) num2str(length(A_cell{i})) ];
    elseif i > 26*3
        range = [ 'C' alphabet(i-26*3) num2str(1) ':C' alphabet(i-26*3) num2str(length(A_cell{i})) ];
    elseif i > 26*2
        range = [ 'B' alphabet(i-26*2) num2str(1) ':B' alphabet(i-26*2) num2str(length(A_cell{i})) ];
    elseif i > 26
        range = [ 'A' alphabet(i-26) num2str(1) ':A' alphabet(i-26) num2str(length(A_cell{i})) ];
    else
        range = [ alphabet(i) num2str(1) ':' alphabet(i) num2str(length(A_cell{i})) ]; 
    end
        

    writematrix(A_cell{i}, filename, 'Range', range );
    % 'Range' specifies the column (A, B, C, etc.) where data is written
    % 'WriteMode' is set to 'append' to add columns for each array
end
