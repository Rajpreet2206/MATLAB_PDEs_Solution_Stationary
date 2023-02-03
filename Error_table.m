function T = into_table(errors, red_factors, N, heading, title)
    headers = [heading , string(N)];
    T = array2table(zeros(0,length(N)+1));
    T.Properties.VariableNames = headers;
    row1 = array2table(["error", errors],'VariableNames',headers);
    row2 = array2table(["error red.", red_factors],'VariableNames',headers);
    T = [T; row1; row2];
    disp(title);
    disp(T);