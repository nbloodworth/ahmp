'''
ahmp.modules.utilities

ahmp module containing esstential application specific utilities.
'''

class Notify:
    '''
    Description:
    A class for formatting output to the command line. Concatenates
    input strings with escape sequences (ANSI formatting codes) 
    to customize output appearance. Helpful for getting the user's
    attention or globally altering message leves.
    '''
    # Define some default values that can be called from any instance for quick custom formatting
    colors={
            "red":'\033[91m',
            "yellow":'\033[93m',
            "green":'\033[92m',
            "blue":'\033[94m',
            "cyan":'\033[96',
            "magenta":'\033[95m',
            "black":'\033[90m',
            "white":'\033[97m'
            }
    formats={
            "END":'\033[0m',
            "bold":'\033[1m',
            "italic":'\033[23m',
            "underline":'\03324m',
            "none":""
            }
    highlights={
            "blue":"\033[48;5;27m",
            "red":"\033[48;5;196m",
            "none":""
            }
    
    _message_level_keys={
            "all":4,
            "notifications":3,
            "warnings":2,
            "errors":1,
            "none":0
            }
    
    _message_format_keys={
            "msg_color":"white",
            "notice_color":"green",
            "warn_color":"yellow",
            "error_color":"red",
            "msg_format":"none",
            "notice_format":"bold",
            "warn_format":"bold",
            "error_format":"bold",
            "msg_highlight":"none",
            "notice_highlight":"none",
            "warn_highlight":"none",
            "error_highlight":"none"
            }

    def __init__(self, toprint="all", **kwargs): 
        '''
        Description:
        Constructor for Notify. Allows specification of message level,
        and custom colors if desired.
        
        Arguments:
        >toprint: specify the message level urgency
            Message levels:
                all         Print all messages
                warnings    Print warnings only
                error       Print errors only
                none        Print no messages
        >**kwargs: other keyword arguments include the components of
        _message_format_key. These are iterated through and set if
        after testing against preset values.
        '''
        # Validate user inputs to the constructor, setting defaults
        # if the user-provided values fall outside allowed parameters
        if toprint not in self._message_level_keys.keys():
            self.toprint=self._message_level_keys["all"]
        else:
            self.toprint=self._message_level_keys[toprint]
        
        # Do the same for the additional keyword arguments
        for k,i in kwargs.items():
            if k in self._message_format_keys.keys():
                if "color" in k and i in self.colors.keys():
                    self._message_format_keys[k]=i
                elif "format" in k and i in self.formats.keys():
                    self._message_format_keys[k]=i
                elif "highlight" in k and i in self.highlights.keys():
                    self._message_format_keys[k]=i
    
    def warning(self,msg,format_only=False,endline="\n"):
        '''
        Description:
        Default method for pre-formatted warnings.

        Arguments:
        >msg: [str] The message to print
        >format_only: [bool][optional] Ignore the value of self.toprint
        and instead return the formatted message as an output
        >endline: [str][optional] argument for print(end=endline)
        '''
        fmsg='{}{}[WARNING]: {}{}{}'.format(
                self.colors[self._message_format_keys["warn_color"]],
                self.highlights[self._message_format_keys["warn_highlight"]],
                self.formats[self._message_format_keys["warn_format"]],
                msg,
                self.formats["END"])
        if format_only:
            return fmsg
        
        if self.toprint>=self._message_level_keys["warnings"]:
            print(fmsg,end=endline)

        return

    def error(self,msg,format_only=False,endline="\n"):
        '''
        Description:
        Default method for pre-formatted errors.

        Arguments:
        >msg: The message to print (string)
        >format_only: [bool][optional] Ignore the value of self.toprint
        and instead return the formatted message as an output
        >endline: [str][optional] argument for print(end=endline)
        '''
        fmsg='{}{}[ERROR]: {}{}{}'.format(
                self.colors[self._message_format_keys["error_color"]],
                self.highlights[self._message_format_keys["error_highlight"]],
                self.formats[self._message_format_keys["error_format"]],
                msg,
                self.formats["END"])
        if format_only:
            return fmsg
        
        if self.toprint>=self._message_level_keys["errors"]:
            print(fmsg,end=endline)

        return
 
    def notice(self,msg,format_only=False, endline="\n"):
        '''
        Description:
        Default method for pre-formatted notices.

        Arguments:
        >msg: The message to print (string)
        >format_only: [bool][optional] Ignore the value of self.toprint
        and instead return the formatted message as an output
        >endline: [str][optional] argument for print(end=endline)
        '''
        fmsg='{}{}{}{}{}'.format(
                self.colors[self._message_format_keys["notice_color"]],
                self.highlights[self._message_format_keys["notice_highlight"]],
                self.formats[self._message_format_keys["notice_format"]],
                msg,
                self.formats["END"])
        if format_only:
            return fmsg
        
        if self.toprint>=self._message_level_keys["notifications"]:
            print(fmsg,end=endline)

        return

 
    def message(self,msg,format_only=False,endline="\n"):
        '''
        Description:
        Print a custom-formatted message.

        Arguments:
        >msg: The message to print (string)
        >format_only: [bool][optional] Ignore the value of self.toprint
        and instead return the formatted message as an output
        >endline: [str][optional] argument for print(end=endline)
        '''
        fmsg='{}{}{}{}{}'.format(
                self.colors[self._message_format_keys["msg_color"]],
                self.highlights[self._message_format_keys["msg_highlight"]],
                self.formats[self._message_format_keys["msg_format"]],
                msg,
                self.formats["END"])
        if format_only:
            return fmsg
        
        if self.toprint>=self._message_level_keys["all"]:
            print(fmsg,end=endline)

        return