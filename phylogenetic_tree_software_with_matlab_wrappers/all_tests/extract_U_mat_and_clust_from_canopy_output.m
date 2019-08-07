function [U, clust]	 = extract_U_mat_and_clust_from_canopy_output(canopy_output)

	canopy_U_matrix = canopy_output{1};
	
	%fprintf('The M_canopy_output{1} matrix:\n');
	disp(canopy_U_matrix);

	num_rows_canopy_U_matrix = length(canopy_U_matrix);
	%fprintf('The number of rows in M_canopy_output{1} matrix = %d\n', num_rows_canopy_U_matrix);

	% Container matrix for new columns being generated
	new_column_silo = [];
	
	% Container matrices for upper diagonal form and LU decomposition
	row_ID_num_1 = [];
	col_ID_num_1 = [];
	col_sorted_canopy_U_matrix = [];
	col_row_sorted_canopy_U_matrix = [];
	row_old_new_position_index = zeros(num_rows_canopy_U_matrix,3);
	upper_diag_canopy_U_matrix = [];

	for i = 1:num_rows_canopy_U_matrix
		U_matrix_row = canopy_U_matrix(i,:);
		%fprintf('The M_canopy_output{1} matrix row being read:\n');
		%disp(U_matrix_row);
		
		if(i ~= 1)
			%disp(sum(U_matrix_row));
			if(sum(U_matrix_row) > 1)
				%new column needs to be created by AND operation
				%finding column positions of the 1s in this row
				num_cols_U_matrix_row = length(U_matrix_row);
				
				%finding column positions of the 1s in this row and storing col
				%positions in array
				col_posn_1_U_matrix_row = find(U_matrix_row == 1);
				num_of_1_U_matrix_row = length(col_posn_1_U_matrix_row);
				%fprintf('The number of 1 in U_matrix_row = %d\n', num_of_1_U_matrix_row);
				%disp(col_posn_1_U_matrix_row);
				
				%Creating a column matrix of 1s and storing in a buffer so that
				%any AND with this, for 1st iteration, will preserve that
				%column
				new_column_generated = ones(num_rows_canopy_U_matrix, 1);
				
				for j = 1:num_of_1_U_matrix_row
					current_column = col_posn_1_U_matrix_row(j);
					new_column_generated = canopy_U_matrix(:,current_column) & new_column_generated;
					%fprintf('new column AND at row: %d and current origin column: %d:\n', i, current_column);
					%disp(new_column_generated);
				end
				%Storing the final generated column, for this row, in a silo
				new_column_silo = [new_column_silo new_column_generated];
				
			end
		end
	end
	
	%fprintf('new columns that can be added to canopy_U_matrix:\n');
	%disp(new_column_silo);
	
	%Concatenating the new columns with the original canopy_U_matrix
	canopy_U_matrix = [canopy_U_matrix new_column_silo];
	%fprintf('The new canopy_U_matrix:\n');
	%disp(canopy_U_matrix);
	%Checking to see if any of these columns have been seen before in the old canopy_U_matrix
	%and deleting silo columns that are similar
	canopy_U_matrix_unique_cols = unique(canopy_U_matrix', 'row', 'stable'); % we are assuming that, when there are repetitions, we delete the copies with higher indices and keep the copies with smaller indices
	canopy_U_matrix = canopy_U_matrix_unique_cols';
	%fprintf('~~~~~~The new canopy_U_matrix after deleting duplicate columns:\n');
	%disp(canopy_U_matrix);
	
	%Calculating the determinant of the canopy_U_matrix
	det_canopy_U_matrix = det(canopy_U_matrix); % if the matrix is not square at this point we get an error
	%fprintf('The determinant of canopy_U_matrix = %d\n', det_canopy_U_matrix);
	
	%If the determinant is 0 then the matrix is not invertible
	if(det_canopy_U_matrix == 0)
		fprintf('The canopy_U_matrix is not invertible!\n');
	else
		fprintf('The canopy_U_matrix is invertible!\n');
		
		num_rows_canopy_U_matrix = length(canopy_U_matrix);
		%fprintf('The number of rows in canopy_U_matrix after column addition = %d\n', num_rows_canopy_U_matrix);
		
		num_cols_canopy_U_matrix = size(canopy_U_matrix,2);
		%fprintf('The number of columns in canopy_U_matrix after column addition = %d\n', num_cols_canopy_U_matrix);
		
		%Transforming canopy_U_matrix into upper diagonal form:
		%for loop to start counting the number of 1s in each column and create
		%a key, value pair of column IDs and number of 1s
		for j_col = 1:num_cols_canopy_U_matrix
			U_matrix_col = canopy_U_matrix(:,j_col);
			%fprintf('The canopy_U_matrix col being read:\n');
			%disp(U_matrix_col);
			
			num_of_1_U_matrix_col = sum(U_matrix_col);
%  			fprintf('The number of 1s in U_matrix_col = %d\n', num_of_1_U_matrix_col);
			
			col_ID_num_1_store = [j_col num_of_1_U_matrix_col];
			col_ID_num_1 = [col_ID_num_1;col_ID_num_1_store];
		end
		
% 		fprintf('The col_ID_num_1 matrix with column IDs and number of 1s:\n');
% 		disp(col_ID_num_1);
		
		%Sorting the col_ID_num_1 matrix according to the number of 1s in 2nd
		%column
		[~,idx] = sort(col_ID_num_1(:,2), 'ascend') %sort the 2nd column, get sort indices
		col_ID_num_1_sorted = col_ID_num_1(idx,:) %sort the whole matrix
% 		fprintf('The sorted col_ID_num_1 matrix according to number of 1s:\n');
% 		disp(col_ID_num_1_sorted);
				
		%for loop to start swapping columns with highest number of 1s in
		%ascending order, right to left
		for j_col = 1:num_cols_canopy_U_matrix
			col_pos_sorted = col_ID_num_1_sorted(j_col,1);
			U_matrix_col = canopy_U_matrix(:,col_pos_sorted);
			
			%Adding columns according to sorted order, descending total 1s
			col_sorted_canopy_U_matrix = [col_sorted_canopy_U_matrix U_matrix_col];
		end
% 		fprintf('The column sorted (descending) canopy_U_matrix matrix according to number of 1s:\n');
% 		disp(col_sorted_canopy_U_matrix);
		
		%for loop to start counting the number of 1s in each row and create
		%a key, value pair of row IDs and number of 1s
		for i_row = 1:num_rows_canopy_U_matrix
			U_matrix_row = col_sorted_canopy_U_matrix(i_row,:);
% 			fprintf('The canopy_U_matrix row being read:\n');
% 			disp(U_matrix_row);
			
			num_of_1_U_matrix_row = sum(U_matrix_row);
% 			fprintf('The number of 1s in U_matrix_row = %d\n', num_of_1_U_matrix_row);
			
			row_ID_num_1_store = [i_row num_of_1_U_matrix_row];
			row_ID_num_1 = [row_ID_num_1;row_ID_num_1_store];
		end
		
% 		fprintf('The row_ID_num_1 matrix with row IDs and number of 1s:\n');
% 		disp(row_ID_num_1);
		
		%Sorting the row_ID_num_1 matrix according to the number of 1s in 2nd
		%column
		[~,idx] = sort(row_ID_num_1(:,2), 'descend') %sort the 2nd column, get sort indices
		row_ID_num_1_sorted = row_ID_num_1(idx,:) %sort the whole matrix
% 		fprintf('The sorted row_ID_num_1 matrix according to number of 1s:\n');
% 		disp(row_ID_num_1_sorted);
		
		%for loop to start swapping rows with highest number of 1s in
		%descending order, top to bottom
		for i_row = 1:num_rows_canopy_U_matrix
			row_pos_sorted = row_ID_num_1_sorted(i_row,1);
			U_matrix_row = col_sorted_canopy_U_matrix(row_pos_sorted,:);
			
			%Adding columns according to sorted order, descending total 1s
			col_row_sorted_canopy_U_matrix = [col_row_sorted_canopy_U_matrix;U_matrix_row];
			
			%Also adding the row positions from row_ID_num_1_sorted to
			%row_old_new_position_index array
			row_old_new_position_index(i_row,1) = row_ID_num_1_sorted(i_row,1);
			row_old_new_position_index(i_row,2) = i_row;
			
		end
% 		fprintf('The row sorted (descending) col_sorted_canopy_U_matrix matrix according to number of 1s:\n');
% 		disp(col_row_sorted_canopy_U_matrix);
% 		fprintf('row_old_new_position_index after row descending arrangement of canopy_U_matrix:\n');
% 		disp(row_old_new_position_index);
		
		matrix_temp_storage = col_row_sorted_canopy_U_matrix;
		% Scan the row_sorted_canopy_U_matrix column by column, and move rows
		% up, so that the i-th column does not have ones below the i-th row
		for j_col = 1:num_cols_canopy_U_matrix
			U_matrix_col = col_row_sorted_canopy_U_matrix(:,j_col);
% 			fprintf('The col_row_sorted_canopy_U_matrix col being read:\n');
% 			disp(U_matrix_col);
			
			%finding row positions of the 1s in this column and storing row
			%positions in array
			row_posn_1_U_matrix_col = find(U_matrix_col == 1);
% 			fprintf('The row_posn_1_U_matrix_col matrix with row positions of 1s:\n');
% 			disp(row_posn_1_U_matrix_col);
			
			%The j_col th column ID is stored at column_position, limit for row
			column_position = j_col;
			
			num_of_1_U_matrix_col = length(row_posn_1_U_matrix_col);
% 			fprintf('The number of 1 in U_matrix_col %d = %d\n', column_position, num_of_1_U_matrix_col);
			
			for i = 1:num_of_1_U_matrix_col
				row_position_1 = row_posn_1_U_matrix_col(i,1);
				if(row_position_1 > column_position)
% 					fprintf('The row_position_1 %d has to be brought up to (best effort) %d\n', row_position_1, column_position);
					%temp = find(row_old_new_position_index(:,2) == row_position_1);
					%row_old_new_position_index(temp, 3) = column_position;
					
					row_store = col_row_sorted_canopy_U_matrix(row_position_1,:);
					r = column_position;
					%Moving the row up to the ith column ID number row, and
					%pushing down the following rows by 1
					b = [col_row_sorted_canopy_U_matrix([1:r-1],:);row_store;col_row_sorted_canopy_U_matrix([r:end],:)];
					%Since all rows were pushed down by the insertion above
					row_delete = row_position_1 + 1;
					b(row_delete,:) = [];
					col_row_sorted_canopy_U_matrix = b;
% 					fprintf('~~~in loop col_row_sorted_canopy_U_matrix, in upper diagonal form:~~~\n');
% 					disp(col_row_sorted_canopy_U_matrix);
					
				elseif(row_position_1 <= column_position)
% 					fprintf('The row_position_1 %d stays at %d\n', row_position_1, row_position_1);
				end
			end
		end
% 		fprintf('The sorted col_row_sorted_canopy_U_matrix, in upper diagonal form:\n');
% 		disp(col_row_sorted_canopy_U_matrix);
% 		fprintf('---row_old_new_position_index of canopy_U_matrix after final transform:\n');
% 		disp(row_old_new_position_index);
		
		%To fill up the 3rd column of row_old_new_position_index in all zero
		%positions
		for i_row = 1:num_rows_canopy_U_matrix
			second_col_row_ID = row_old_new_position_index(i_row,2);
			matrix_temp_storage_row = matrix_temp_storage(second_col_row_ID,:);

			%matching the matrix_temp_storage_row with a row in the col_row_sorted_canopy_U_matrix
			for j = 1:num_rows_canopy_U_matrix
				if(all(matrix_temp_storage_row == col_row_sorted_canopy_U_matrix(j,:)) == 1)
					row_old_new_position_index(i_row,3) = j;
					break
				end
			end
		end
% 		fprintf('---row_old_new_position_index of canopy_U_matrix after final transform and 3rd column:\n');
% 		disp(row_old_new_position_index);
	end
	
	%1st column of row_old_new_position_index is the row positions for canopy_U_matrix
	%3rd columns are the corresponding new row positions for
	%col_row_sorted_canopy_U_matrix which will be U2
	%Reordering cluster IDs according to this mapping to get clust2
	%canopy_clonalMut_cluster = canopy_output{6};
	
	%
	final_clust = [];
	for i = 1:num_rows_canopy_U_matrix
		tmp = canopy_output{3}{  row_old_new_position_index(i,1)  };
		for j = tmp
			final_clust = [final_clust ; [ j  , row_old_new_position_index(i,3)  ]];
		end
	end
	
	clust = final_clust;
	
	U = col_row_sorted_canopy_U_matrix;
	
end