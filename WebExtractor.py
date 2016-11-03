import time
import glob
import urllib

from ObjectBase import *
from FileHandler import *


FULL_WRITER = 'full'
VI_WRITER = 'vi'


class Journal(SimpleObjectBase):
	def __init__(self, name, date, volume, pages, impactFactor = None, issnNb = None, aliasList = None):
		init(self, locals())
	
	@staticmethod
	def getJournalDateVolumeAndPagesFromStr(journalInfoStr):
		matchObj = re.match('^([\w\s\-\(\)\[\]]+)\. ([\w\d\s\-]+)(;([\w\-\(\)\s]+):(.*)\.)?', journalInfoStr)
		if matchObj is None:
			if 'ahead of' in journalInfoStr:
				return journalInfoStr.split(' [')[0], Utilities.getShortTimeString()[:4], '', ''
			print 'Journal info [%s] not handled' % journalInfoStr
			raise NotImplementedError
		journalName = matchObj.group(1)
		if '.' in journalName:
			print 'Journal [%s] is weird... [%s]' % (journalName, journalInfoStr)
			raise NotImplementedError
		date = matchObj.group(2)
		volume = pages = None
		if len(matchObj.groups()) > 3:
			volume = matchObj.group(4)
			pages = matchObj.group(5)
		return journalName, date, volume, pages
	
	def getJournalInfo(self):
		journalInfoStr = '%s. %s' % (self.name, self.date)
		if self.volume and self.pages:
			journalInfoStr += '; %s:%s' % (self.volume, self.pages)
		return journalInfoStr
		

class Article(SimpleObjectBase):
	PMID = 'PMID'
	LEVEL1DIET_ID = 'LV1ID'
	
	def __init__(self, idList, title, abstract, authors, institute, journal, keywords = None, linkList = None,
				 publishedDate = None, key = None):
		if type(idList) != types.ListType:
			idList = [idList]
		init(self, locals())
	
	def getImpactFactor(self):
		return ['%0.2f' % impactFactor for impactFactor in self.journal.impactFactor]
	
	def getPubmedId(self):
		idDict = dict([[idName, id] for idName, id in self.idList])
		return idDict.get(self.PMID)
	
	def getJournalInfo(self):
		return self.journal.getJournalInfo()
	
	def getId(self):
		return ', '.join([':'.join(id) for id in self.idList])
		

class WebExtractor(object):
	_endPage = None
	
	def __init__(self, endPage = None, nbRetries = 5):
		self._nbRetries = nbRetries
		if endPage is None and self._endPage is None:
			endPage = '</html>'
		self._endPage = endPage.lower()
	
	def _getStrIncludedInTag(self, text, tag, endTag = None):
		startIdx = text.find(tag)
		if endTag is None:
			endTag = tag[0] + '/' + tag[1:]
		self._endIdx = text[startIdx+len(tag):].find(endTag)
		#print 1, [tag, startIdx, endTag, self._endIdx, text[startIdx+len(tag): startIdx+100]]
		if startIdx == -1 or self._endIdx == -1:
			return
		self._endIdx += startIdx + len(tag)
		#print 2, [tag, startIdx, endTag, self._endIdx, text[startIdx+len(tag): startIdx+100]]
		return text[startIdx+len(tag): self._endIdx]
	
	def _getStrListIncludedInTag(self, text, tag, endTag = None):
		strList = []
		while True:
			extractedStr = self._getStrIncludedInTag(text, tag, endTag)
			if extractedStr is None:
				break
			strList.append(extractedStr)
			text = text[self._endIdx+len(endTag):]
		return strList
	
	def _saveUrlToFile(self, url, fileName):
		content = url
		if not '\n' in url:
			content = WebExtractor()._getPageContentForUrl(url)
		fh = open(fileName, 'w')
		fh.write(content)
		fh.close()
	
	def _getPageContentForUrl(self, url):
		import urllib
		for i in xrange(self._nbRetries):
			fh = urllib.urlopen(url)
			content = fh.read()
			fh.close()
			endKeyword = content.strip()[-len(self._endPage):].lower()
			if self._endPage == '' or endKeyword == self._endPage:
				return content
			print 'Retrying [%s]: endPage should be "%s" but "%s" found' % (url, self._endPage, [endKeyword])
	
class ArticleWriterBase:
	def _getArticleListFromArticleDict(self, articleDict):
		articleList = []
		for keyword, currentArticleList in articleDict.items():
			articleList += currentArticleList
		return articleList
	
	def _write(self, fileName, articleList):
		raise NotImplemented
	
	def write(self, fileName, articleDict):
		articleList = self._getArticleListFromArticleDict(articleDict)
		self._write(fileName, articleList)
	
	
class EasilyArticleReadableInViWriter(ArticleWriterBase):
	def _write(self, fileName, articleList):
		articleList.sort(lambda article1, article2: cmp((article1.publishedDate, article1.title), (article2.publishedDate, article2.title)))
		fh = open(fileName, 'w')
		fh.write('%d articles extracted\n\n' % len(articleList))
		for article in articleList:
			fh.write('keyword: %s\npubmedId: %d\npublishedDate: %s\nauthors: %s\nTITLE: %s\njournal: %s\nimpact factor: %s\n\n%s\n\n\n' % \
					 (article.key, int(article.getPubmedId()), article.publishedDate, article.authors,
					  article.title, article.getJournalInfo(), article.getImpactFactor(), article.abstract))
		fh.close()
		

class FullArticleWriter(ArticleWriterBase):
	def _write(self, fileName, articleList):
		articleList.sort(lambda article1, article2: cmp(article2.journal.impactFactor, article1.journal.impactFactor))
		fh = CsvFileWriter(fileName)
		fh.write(['Keyword', 'Title', 'Autors', 'Journal', 'Impact factor', 'Institute', 'Links', 'Abstract'])
		for article in articleList:
			fh.write([article.key, article.title, article.authors, article.getJournalInfo(),
					  article.getImpactFactor(), article.institute, ','.join(article.linkList), article.abstract])

			
class ArticleExtractorBase(object):
	sleepTimeInSec = 2
	endPage = None
	
	def __init__(self, saveProcessedIdList = False):
		self._webExtractor = WebExtractor(endPage = self.endPage)
		self._processedIdList = []
		self._saveProcessedIdList = saveProcessedIdList
		self._articleNb = 0
	
	def _getStrIncludedInTag(self, content, startTag, endTag = None):
		resStr = self._webExtractor._getStrIncludedInTag(content, startTag, endTag)
		self._endIdx = self._webExtractor._endIdx
		return resStr
		
	def _getPageContentForUrl(self, url):
		return self._webExtractor._getPageContentForUrl(url)
		
	def _getArticleFromPageContent(self, pageContent, url):
		raise NotImplementedError
		
	def _getArticleByUrl(self, url):
		try:
			pageContent = self._getPageContentForUrl(url)
			return self._getArticleFromPageContent(pageContent, url)
		except:
			print 'Failed to parse url [%s]' % url
			raise
	
	def _getUrlListForKeyword(self, keyword):
		raise NotImplementedError
	
	def _getIdFromUrl(self, url):
		return url
	
	def getArticleListByKeyword(self, keyword):
		urlList = self._getUrlListForKeyword(keyword)
		print '%d article(s) extracted for keyword [%s]' % (len(urlList), keyword)
		articleList = []
		for url in urlList:
			currentId = self._getIdFromUrl(url)
			if currentId in self._processedIdList:
				print 'Id [%s] already processed' % str(currentId)
				continue
			self._processedIdList.append(currentId)
			article = self._getArticleByUrl(url)
			if article:
				article.key = keyword
				articleList.append(article)
			time.sleep(self.sleepTimeInSec)
		return articleList
	
	def _setAlreadyProcessedIdList(self):
		return
	
	def getArticleDictForKeywordList(self, keywordList, dirName = None):
		self._dirName = dirName
		self._processedIdList = []
		self._setAlreadyProcessedIdList()
		keywordToArticleDict = {}
		for keyword in keywordList:
			print '-' * 100
			print 'Processing keyword [%s]' % keyword
			articleList = self.getArticleListByKeyword(keyword)
			self._articleNb += len(articleList)
			keywordToArticleDict[keyword] = articleList
		return keywordToArticleDict
	
	def _getWriterClassFromShortName(self, writerShortName):
		writerClass = articleWriterDict.get(writerShortName)
		if writerClass is None:
			print 'Writer class should be in %s' % articleWriterDict.keys()
			raise NotImplementedError
		return writerClass
	
	def _getFileNameFromDirAndFileName(self, dirName, fileName):
		if dirName:
			fileName = os.path.join(dirName, fileName)
		return fileName
	
	def process(self, keywordList, fileName = 'abstracts', dirName = None, writerShortName = FULL_WRITER):
		writerClass = self._getWriterClassFromShortName(writerShortName)
		if dirName and not Utilities.doesDirExist(dirName):
			os.mkdir(dirName)
		self._articleNb = 0
		self.__dateStr = Utilities.getShortTimeString()
		self._fileName = fileName
		articleDict = self.getArticleDictForKeywordList(keywordList, dirName)
		print 'total: %d articles extracted' % self._articleNb
		self.writeDict(self._getFileNameFromDirAndFileName(dirName, self._fileName + '_%s.csv' % self.__dateStr),
					   articleDict, writerClass)
		if self._saveProcessedIdList:
			Utilities.saveCache(self._processedIdList, self._getFileNameFromDirAndFileName(dirName, fileName + '_%s.pyDump' % self.__dateStr))
	
	def writeDict(self, fileName, articleDict, writerClass):
		writerClass().write(fileName, articleDict)
	
articleWriterDict = {'vi': EasilyArticleReadableInViWriter, 'full': FullArticleWriter}
