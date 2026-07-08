function methodName = normalizeMethodName(methodName)
% normalizeMethodName  Validate and normalize the reconstruction method name.

methodName = upper(strtrim(string(methodName)));
if methodName == "WBARL"
    methodName = "WB-ARL";
end

if methodName ~= "RLD" && methodName ~= "WB-ARL"
    error('Unsupported methodName "%s". Use "RLD" or "WB-ARL".', methodName);
end

end
