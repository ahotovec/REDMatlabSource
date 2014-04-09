function saveTemplSearch(D,monthvec,staCode)
save(['C:\Users\geo-user\Documents\KateA\MultCodes\RCStemplsearch\',datestr(monthvec,'yyyy-mm_'),staCode,'.mat'],'D');
end
