''' Team: - Prapat Suriyaphol
          - Masao Yamaguchi
          - Victor Renault

    Research project: Trans-ethnic Genomics Study of Multigenic Disorders
    by Japan/France International Collaboration
	
    Project: Set up a database based system to help the study
			
    Research director: Dr. Fumihiko Matsuda '''


class Color:
	''' defines the basic operations for HTML colors '''
	
	MAX_COLOR = pow(256, 3) - 1
	RED_MASK = 255
	GREEN_MASK = 65280
	BLUE_MASK = 16711680
	HEX1_MASK = 15
	HEX2_MASK = 240
	
	''' nb of colors that gives a good color range '''
	nbColorList = [60, 180]
	
	hexaDict = {10: 'A', 11: 'B', 12: 'C', 13: 'D', 14: 'E', 15: 'F'}
	
	def __getHexadecimalCharFromNumber(self, hexaNb):
		return '%s' % self.hexaDict.get(hexaNb, hexaNb)
	
	def __getHexaDecimalCodeFromNumber(self, colorNumber):
		firstHexValue = colorNumber & self.HEX1_MASK
		secondHexValue = (colorNumber & self.HEX2_MASK) / pow(2, 4)
		return self.__getHexadecimalCharFromNumber(firstHexValue) + self.__getHexadecimalCharFromNumber(secondHexValue)
	
	def __getHtmlColorCode(self, colorNb):
		redValue = colorNb & self.RED_MASK
		greenValue = (colorNb & self.GREEN_MASK) / pow(2, 8)
		blueValue = (colorNb & self.BLUE_MASK) / pow(2, 16)
		return '#' + self.__getHexaDecimalCodeFromNumber(redValue) + self.__getHexaDecimalCodeFromNumber(greenValue) + \
			   self.__getHexaDecimalCodeFromNumber(blueValue)
		
	def __getCurrentColorNb(self, nbColors):
		for currentColorNb in self.nbColorList:
			if currentColorNb < nbColors:
				continue
			return currentColorNb
	
	def getColorList(self, nbColors):
		currentColorNb = self.__getCurrentColorNb(nbColors)
		colorStep = self.MAX_COLOR / (currentColorNb + 1)
		currentColor = colorStep
		colorList = []
		for i in range(nbColors):
			colorList.append(self.__getHtmlColorCode(currentColor))
			currentColor += colorStep
		return colorList
	
	def __getOppositeHexChar(self, hexChar):
		oppositeHexDict = {'F': 0, 'E': 1, 'D': 2, 'C': 3, 'B': 4, 'A': 5, '9': 6, '8': 7, '7': 8, '6': 9, '5': 'A', '4': 'B', '3': 'C', '2': 'D',
						   '1': 'E', '0': 'F'}
		return str(oppositeHexDict[hexChar])
	
	def getOppositeColor(self, htmlColorNbInHexadecimal):
		oppositeColor = '#'
		for i in range(1, 7):
			oppositeColor += self.__getOppositeHexChar(htmlColorNbInHexadecimal[i])
		return oppositeColor


