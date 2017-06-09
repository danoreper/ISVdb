
run = function()
{
    source("./lm/formulaWrapper.R")
    s1 = "~ 1 + x + X + z:y + X:Y + X:Y:Z + (1+ X |RIX) + (1 + X:Y|RIX:Diet) + (Y+1|X)" 
    formulaWrapper$parseCovariateString(s1)

    print(s1)
    ##print(formulaWrapper$removeEffectAndInteractions("X",s1)$modifiedString)
    print(formulaWrapper$onlyKeepEffectAndInteractions(c("1", "X"), s1)$modifiedString)
    ##print(formulaWrapper$onlyKeepEffectAndInteractions(c("X:Y"), s1)$modifiedString)
}
