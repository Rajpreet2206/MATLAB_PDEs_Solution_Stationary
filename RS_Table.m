function T = RS_Table(L, S, N, heading, title)
    headers = [heading , string(N)];
    T = array2table(zeros(0,length(N)+1));
    T.Properties.VariableNames = headers;
    row1 = array2table(["Runtime", L],'VariableNames',headers);
    row2 = array2table(["Storage(count of doubles used)", S],'VariableNames',headers);
    %row3 = array2table(["stability", stability],'VariableNames',headers);
    T = [T; row1; row2];
    disp(title);
    disp(T);