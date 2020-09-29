% function [nnz, colIDs, rowDel] = mat2csr(inMat)
%     [vertices, ~] = size(inMat); 
%     elementsCount = sum(sum(inMat ~= 0));
%     nnzPtr = 1;
%     nnz = zeros(1, elementsCount);
%     colIDs = zeros(1, elementsCount);
%     rowDel = zeros(1, vertices + 1);
%     for i=1:vertices
%         for j=1:vertices
%             if(inMat(i, j) ~= 0)
%                nnz(nnzPtr) = inMat(i, j);
%                colIDs(nnzPtr) = j - 1;
%                nnzPtr = nnzPtr+1;
%             end
%         end
%         rowDel(i+1) = nnzPtr - 1;
%     end
%     
%     csvwrite('csr.csv', [nnz; colIDs]);
%     dlmwrite('csr.csv',rowDel,'-append');
% end

function [nnz, colIDs, rowDel] = mat2csr(inMat)
      [rows,colIDs,nnz] = find(inMat);
      % Sort by row IDs
      A = sortrows([rows, colIDs, nnz])';
      colIDs = A(2,:);
      nnz = ones(1, size(colIDs, 2));
      rowDel = zeros(1, A(1, end) + 1);
      rowDel(1) = 1;
      rowDel(A(1, end) + 1) = size(colIDs, 2) + 1;
      rowDelIdx = 2;
      % Form row delimiters array
      for i=2:length(A(1,:))
          if(A(1, i) ~= A(1, i-1))
              rowDel(rowDelIdx) = i;
              rowDelIdx = rowDelIdx + 1;
          end
      end
      file_id = fopen('csr.csv','wt');
      fprintf(file_id,'%d,', nnz);
      fprintf(file_id, '\n');
      fprintf(file_id,'%d,', colIDs - 1);
      fprintf(file_id, '\n');
      fprintf(file_id,'%d,', rowDel - 1);
      fprintf(file_id, '\n');
      fclose(file_id);
      
      fprintf('GRAPH_DIM %d\n', size(inMat, 1));
      fprintf('nnz %d\n', size(colIDs, 2));
end