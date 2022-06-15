function [ outstring ] = nameParenthesisPreslasher( instring )
%SLW wrote this function to put an escape slash before parenthesis because
%some bash functions require that.
    outstring='';letstore='a';
    for i=1:length(instring)
        let=instring(i);
        if (let=='(' | let==')') && letstore~='\'
            outstring=[outstring, '\'];
        end
        outstring=[outstring,let];
        letstore=let;
    end
    instring=outstring;

end

