no1 = '30'

modelNode1 = getNode(no1+'_right')
modelNode2 = getNode(no1+'_prediction')
modelNode2.SetSpacing(modelNode1.GetSpacing())
modelNode2.SetOrigin(modelNode1.GetOrigin())
